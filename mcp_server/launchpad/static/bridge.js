/**
 * STAR Launchpad ↔ editor bridge client.
 *
 * Connects to a paired VS Code / Cursor extension over localhost WebSocket,
 * authenticates with a shared token stored in localStorage, and exposes a
 * minimal API on an Alpine component instance (launchpadApp). Phase 1 only
 * covers pairing + auth round-trip; submit / streams / cancel stubs exist
 * so UI wiring can land before LM integration in phase 2.
 *
 * Editors persist under key "launchpadBridge.editors" as:
 *   { version: 1, editors: [ { id, label, port, token, host, workspace, addedAt } ] }
 *
 * Only same-origin loopback connections are expected (the extension binds to
 * 127.0.0.1). The token is required as the first message within 2s or the
 * server closes the socket.
 */

(function () {
  const STORE_KEY = "launchpadBridge.editors";
  const STORE_VERSION = 1;
  const RECONNECT_BASE_MS = 1000;
  const RECONNECT_CAP_MS = 30_000;
  const AUTH_SEND_DELAY_MS = 50; // give the socket a tick before sending

  function loadEditors() {
    let raw = null;
    try {
      raw = localStorage.getItem(STORE_KEY);
    } catch {
      return [];
    }
    if (!raw) return [];
    try {
      const parsed = JSON.parse(raw);
      if (!parsed || parsed.version !== STORE_VERSION) return [];
      return Array.isArray(parsed.editors) ? parsed.editors : [];
    } catch {
      return [];
    }
  }

  function saveEditors(editors) {
    try {
      localStorage.setItem(
        STORE_KEY,
        JSON.stringify({ version: STORE_VERSION, editors })
      );
    } catch {
      /* localStorage unavailable (private mode?) — ignore */
    }
  }

  function genEditorId() {
    if (typeof crypto !== "undefined" && crypto.randomUUID) {
      return crypto.randomUUID();
    }
    return "ed-" + Math.random().toString(36).slice(2, 10) + Date.now().toString(36);
  }

  /** Parse the current URL for ?pair=1&port=...&token=... and return params or null. */
  function readPairParams() {
    const sp = new URLSearchParams(window.location.search);
    if (sp.get("pair") !== "1") return null;
    const port = parseInt(sp.get("port") || "", 10);
    const token = sp.get("token") || "";
    if (!Number.isInteger(port) || port <= 0 || !/^[0-9a-f]{32}$/.test(token)) {
      return null;
    }
    return {
      port,
      token,
      label: sp.get("label") || "Editor",
    };
  }

  function stripPairParams() {
    try {
      const url = new URL(window.location.href);
      ["pair", "port", "token", "label"].forEach((k) => url.searchParams.delete(k));
      window.history.replaceState({}, "", url.toString());
    } catch {
      /* ignore */
    }
  }

  /**
   * Bridge client bound to an Alpine component. The component provides the
   * reactive state fields:
   *   bridgeEditors: []
   *   bridgeSelectedId: string | null
   *   bridgeConnection: { status, app, workspace, host, error } | null
   *   bridgeNotice: string
   *   bridgeStreams: Record<id, string>
   */
  class BridgeClient {
    constructor(component) {
      this.c = component;
      this.ws = null;
      this.currentEditor = null;
      this.reconnectAttempt = 0;
      this.reconnectTimer = null;
      this.manualClose = false;
      this.pendingIds = new Set();
    }

    /** Add an editor from a pair URL (dedupes by host+port+token). */
    registerEditor({ port, token, label }) {
      const editors = loadEditors();
      const match = editors.find(
        (e) => e.port === port && e.token === token
      );
      const entry = match || {
        id: genEditorId(),
        label,
        port,
        token,
        host: "",
        workspace: "",
        addedAt: Date.now(),
      };
      if (!match) editors.push(entry);
      saveEditors(editors);
      this._pushToComponent(editors);
      return entry;
    }

    removeEditor(id) {
      const editors = loadEditors().filter((e) => e.id !== id);
      saveEditors(editors);
      if (this.c.bridgeSelectedId === id) {
        this.disconnect();
        this.c.bridgeSelectedId = null;
      }
      this._pushToComponent(editors);
    }

    refresh() {
      this._pushToComponent(loadEditors());
    }

    /** Connect to editor by id. Reconnects with exponential backoff on close. */
    connect(id) {
      const editors = loadEditors();
      const editor = editors.find((e) => e.id === id);
      if (!editor) {
        this.c.bridgeConnection = {
          status: "error",
          error: "Editor not found in local store.",
        };
        return;
      }
      this.disconnect();
      this.currentEditor = editor;
      this.c.bridgeSelectedId = id;
      this.manualClose = false;
      this.reconnectAttempt = 0;
      this._open();
    }

    disconnect() {
      this.manualClose = true;
      if (this.reconnectTimer) {
        clearTimeout(this.reconnectTimer);
        this.reconnectTimer = null;
      }
      if (this.ws) {
        try {
          this.ws.close();
        } catch {
          /* ignore */
        }
        this.ws = null;
      }
      this.currentEditor = null;
      this.c.bridgeConnection = null;
    }

    /** Phase 1 stub; returns the generated id or null on failure. */
    submit(mode, formData, context) {
      if (!this.ws || this.ws.readyState !== WebSocket.OPEN) {
        this.c.bridgeNotice = "Bridge is not connected.";
        return null;
      }
      if (this.c.bridgeConnection?.status !== "connected") {
        this.c.bridgeNotice = "Bridge is not authenticated yet.";
        return null;
      }
      const id = genEditorId();
      this.pendingIds.add(id);
      this.c.bridgeStreams = { ...this.c.bridgeStreams, [id]: "" };
      try {
        this.ws.send(
          JSON.stringify({
            type: "submit",
            id,
            mode,
            formData: formData || {},
            context: context || [],
          })
        );
      } catch (e) {
        this.c.bridgeNotice = "Bridge send failed: " + (e?.message || e);
        this.pendingIds.delete(id);
        return null;
      }
      return id;
    }

    cancel(id) {
      if (!this.ws || this.ws.readyState !== WebSocket.OPEN) return;
      try {
        this.ws.send(JSON.stringify({ type: "cancel", id }));
      } catch {
        /* ignore */
      }
    }

    // --- internals -----------------------------------------------------

    _open() {
      const ed = this.currentEditor;
      if (!ed) return;
      this.c.bridgeConnection = { status: "connecting" };
      let ws;
      try {
        ws = new WebSocket(`ws://127.0.0.1:${ed.port}`);
      } catch (e) {
        this._onFailure(e?.message || String(e));
        return;
      }
      this.ws = ws;

      ws.addEventListener("open", () => {
        setTimeout(() => {
          if (ws.readyState !== WebSocket.OPEN) return;
          try {
            ws.send(JSON.stringify({ type: "auth", token: ed.token }));
            this.c.bridgeConnection = { status: "authenticating" };
          } catch (e) {
            this._onFailure(e?.message || String(e));
          }
        }, AUTH_SEND_DELAY_MS);
      });

      ws.addEventListener("message", (evt) => this._onMessage(evt));
      ws.addEventListener("close", (evt) => this._onClose(evt));
      ws.addEventListener("error", () => {
        // The 'close' handler runs afterwards and does the reconnect work.
        this.c.bridgeNotice = "Bridge socket error.";
      });
    }

    _onMessage(evt) {
      let msg;
      try {
        msg = JSON.parse(evt.data);
      } catch {
        return;
      }
      switch (msg.type) {
        case "auth_ok": {
          this.reconnectAttempt = 0;
          this.c.bridgeConnection = {
            status: "connected",
            app: msg.app,
            workspace: msg.workspace,
            host: msg.host,
            protocol: msg.protocol,
          };
          // Mirror the editor's labels into the stored record.
          const editors = loadEditors();
          const ed = editors.find((e) => e.id === this.currentEditor?.id);
          if (ed) {
            ed.host = msg.host || ed.host;
            ed.workspace = msg.workspace || ed.workspace;
            if (!ed.label || ed.label === "Editor") ed.label = msg.app || ed.label;
            saveEditors(editors);
            this._pushToComponent(editors);
          }
          break;
        }
        case "auth_fail": {
          this.manualClose = true; // don't loop on bad token
          this.c.bridgeConnection = {
            status: "auth_failed",
            error: msg.reason || "auth failed",
          };
          break;
        }
        case "ack": {
          // Acked; nothing to display yet.
          break;
        }
        case "chunk": {
          const prev = this.c.bridgeStreams[msg.id] || "";
          this.c.bridgeStreams = {
            ...this.c.bridgeStreams,
            [msg.id]: prev + (msg.text || ""),
          };
          break;
        }
        case "complete": {
          this.pendingIds.delete(msg.id);
          break;
        }
        case "error": {
          this.pendingIds.delete(msg.id);
          this.c.bridgeNotice = msg.message || "bridge error";
          break;
        }
        case "editor_event": {
          // Phase 1: just record the latest event; no UI yet.
          this.c.bridgeLastEvent = {
            kind: msg.kind,
            payload: msg.payload,
            at: Date.now(),
          };
          break;
        }
        default:
          break;
      }
    }

    _onClose(evt) {
      this.ws = null;
      if (this.manualClose) {
        this.c.bridgeConnection = null;
        return;
      }
      const ok = evt && evt.code === 1000;
      this.c.bridgeConnection = {
        status: ok ? "closed" : "reconnecting",
        error: evt && evt.reason ? evt.reason : "",
      };
      this._scheduleReconnect();
    }

    _onFailure(msg) {
      this.ws = null;
      this.c.bridgeConnection = { status: "error", error: msg };
      this._scheduleReconnect();
    }

    _scheduleReconnect() {
      if (this.manualClose || !this.currentEditor) return;
      const delay = Math.min(
        RECONNECT_BASE_MS * Math.pow(2, this.reconnectAttempt),
        RECONNECT_CAP_MS
      );
      this.reconnectAttempt += 1;
      this.reconnectTimer = setTimeout(() => {
        this.reconnectTimer = null;
        this._open();
      }, delay);
    }

    _pushToComponent(editors) {
      this.c.bridgeEditors = editors.slice();
    }
  }

  function installBridge(component) {
    const client = new BridgeClient(component);
    component.$bridge = client;
    component.bridgeEditors = loadEditors();

    // Consume /pair query params if present.
    const pair = readPairParams();
    if (pair) {
      const entry = client.registerEditor(pair);
      component.bridgeNotice =
        `Paired with "${pair.label}" on port ${pair.port}. Select it below to connect.`;
      component.bridgeSelectedId = entry.id;
      stripPairParams();
      // Auto-connect after Alpine has rendered.
      setTimeout(() => client.connect(entry.id), 150);
    }
  }

  window.installLaunchpadBridge = installBridge;
  window.loadLaunchpadEditors = loadEditors;
})();
