import * as net from "net";
import * as os from "os";
import * as vscode from "vscode";
import { WebSocketServer, WebSocket } from "ws";
import {
  ClientMessage,
  PROTOCOL_VERSION,
  ServerMessage,
} from "./protocol";
import { tokenMatches } from "./token";

const AUTH_TIMEOUT_MS = 2000;
const LOOPBACK_HOST = "127.0.0.1";

export type BridgeState =
  | { kind: "stopped" }
  | { kind: "listening"; port: number; clients: number }
  | { kind: "error"; message: string };

export interface BridgeEvents {
  onStateChanged: vscode.Event<BridgeState>;
  onClientCountChanged: vscode.Event<number>;
}

/** Wraps a ws WebSocketServer bound to 127.0.0.1 with token auth. */
export class BridgeServer implements vscode.Disposable {
  private wss: WebSocketServer | null = null;
  private port: number | null = null;
  private authed = new Set<WebSocket>();
  private readonly stateEmitter = new vscode.EventEmitter<BridgeState>();
  private readonly clientEmitter = new vscode.EventEmitter<number>();

  readonly onStateChanged = this.stateEmitter.event;
  readonly onClientCountChanged = this.clientEmitter.event;

  constructor(
    private readonly getToken: () => string,
    private readonly logger: vscode.OutputChannel
  ) {}

  get currentPort(): number | null {
    return this.port;
  }

  get clientCount(): number {
    return this.authed.size;
  }

  /** Pick an open port in [low, high], bind, and start listening. */
  async start(range: [number, number]): Promise<void> {
    if (this.wss) {
      return;
    }
    const port = await this.pickPort(range);
    this.logger.appendLine(`Binding bridge WebSocket on 127.0.0.1:${port}`);
    const wss = new WebSocketServer({ host: LOOPBACK_HOST, port });
    this.wss = wss;
    this.port = port;

    wss.on("connection", (socket) => this.handleConnection(socket));
    wss.on("error", (err) => {
      this.logger.appendLine(`WebSocketServer error: ${err.message}`);
      this.stateEmitter.fire({ kind: "error", message: err.message });
    });

    this.stateEmitter.fire({ kind: "listening", port, clients: 0 });
  }

  /** Send a message to every authed client (used by editor_event fanout). */
  broadcast(msg: ServerMessage): void {
    const payload = JSON.stringify(msg);
    for (const c of this.authed) {
      if (c.readyState === WebSocket.OPEN) {
        c.send(payload);
      }
    }
  }

  /** Force-close every open socket (e.g. after resetToken). */
  invalidateAllSessions(reason: string): void {
    if (!this.wss) return;
    for (const client of this.wss.clients) {
      try {
        this.safeSend(client, { type: "auth_fail", reason });
        client.close(1008, reason);
      } catch {
        /* ignore */
      }
    }
    this.authed.clear();
    this.clientEmitter.fire(0);
  }

  dispose(): void {
    if (this.wss) {
      this.wss.close();
      this.wss = null;
    }
    this.port = null;
    this.authed.clear();
    this.stateEmitter.fire({ kind: "stopped" });
    this.stateEmitter.dispose();
    this.clientEmitter.dispose();
  }

  // --- internals ---------------------------------------------------------

  private handleConnection(socket: WebSocket): void {
    this.logger.appendLine("Bridge client connected; awaiting auth.");
    const authTimer = setTimeout(() => {
      if (!this.authed.has(socket)) {
        this.safeSend(socket, { type: "auth_fail", reason: "auth timeout" });
        socket.close(1008, "auth timeout");
      }
    }, AUTH_TIMEOUT_MS);

    socket.on("message", (raw) => {
      let msg: ClientMessage;
      try {
        msg = JSON.parse(raw.toString("utf8")) as ClientMessage;
      } catch {
        this.safeSend(socket, {
          type: "error",
          id: "",
          message: "Invalid JSON payload",
        });
        return;
      }
      if (!this.authed.has(socket)) {
        this.handleUnauthenticated(socket, msg, authTimer);
        return;
      }
      this.handleAuthenticated(socket, msg);
    });

    socket.on("close", () => {
      clearTimeout(authTimer);
      if (this.authed.delete(socket)) {
        this.clientEmitter.fire(this.authed.size);
      }
    });

    socket.on("error", (err) => {
      this.logger.appendLine(`Bridge socket error: ${err.message}`);
    });
  }

  private handleUnauthenticated(
    socket: WebSocket,
    msg: ClientMessage,
    authTimer: NodeJS.Timeout
  ): void {
    if (msg.type !== "auth") {
      this.safeSend(socket, {
        type: "auth_fail",
        reason: "expected auth as first message",
      });
      socket.close(1008, "auth expected");
      return;
    }
    if (!tokenMatches(this.getToken(), msg.token)) {
      this.safeSend(socket, { type: "auth_fail", reason: "token mismatch" });
      socket.close(1008, "token mismatch");
      return;
    }
    clearTimeout(authTimer);
    this.authed.add(socket);
    this.clientEmitter.fire(this.authed.size);
    this.safeSend(socket, {
      type: "auth_ok",
      app: vscode.env.appName ?? "unknown",
      workspace: vscode.workspace.name ?? "",
      host: os.hostname(),
      protocol: PROTOCOL_VERSION,
    });
    this.logger.appendLine("Bridge client authenticated.");
  }

  private handleAuthenticated(socket: WebSocket, msg: ClientMessage): void {
    switch (msg.type) {
      case "submit": {
        // Phase 1: LM wiring lands in phase 2. Ack + stub response so the SPA
        // can verify the full round-trip today.
        this.safeSend(socket, { type: "ack", id: msg.id });
        const stubText =
          `[bridge phase 1] received ${msg.mode} request with ` +
          `${Object.keys(msg.formData || {}).length} form field(s) and ` +
          `${(msg.context || []).length} context message(s). ` +
          `LM streaming lands in phase 2.`;
        this.safeSend(socket, { type: "chunk", id: msg.id, text: stubText });
        this.safeSend(socket, { type: "complete", id: msg.id });
        break;
      }
      case "followup":
      case "cancel": {
        this.safeSend(socket, {
          type: "error",
          id: msg.id,
          message: `'${msg.type}' lands in phase 2`,
        });
        break;
      }
      case "auth": {
        // Already authed; ignore repeat auth.
        break;
      }
      default: {
        this.safeSend(socket, {
          type: "error",
          id: "",
          message: `Unknown message type: ${(msg as { type: string }).type}`,
        });
      }
    }
  }

  private safeSend(socket: WebSocket, msg: ServerMessage): void {
    if (socket.readyState === WebSocket.OPEN) {
      socket.send(JSON.stringify(msg));
    }
  }

  private async pickPort(range: [number, number]): Promise<number> {
    const [low, high] = range;
    if (!Number.isInteger(low) || !Number.isInteger(high) || low > high) {
      throw new Error(`Invalid port range [${low}, ${high}]`);
    }
    for (let p = low; p <= high; p++) {
      if (await this.isPortFree(p)) {
        return p;
      }
    }
    throw new Error(`No free TCP port in [${low}, ${high}] on 127.0.0.1`);
  }

  private isPortFree(port: number): Promise<boolean> {
    return new Promise((resolve) => {
      const tester = net.createServer();
      tester.once("error", () => resolve(false));
      tester.once("listening", () => {
        tester.close(() => resolve(true));
      });
      tester.listen(port, LOOPBACK_HOST);
    });
  }
}
