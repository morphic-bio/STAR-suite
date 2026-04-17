import * as vscode from "vscode";
import { BridgeServer, BridgeState } from "./bridge-server";
import { getOrCreateToken, rotateToken } from "./token";

let server: BridgeServer | null = null;
let statusBar: vscode.StatusBarItem | null = null;
let logger: vscode.OutputChannel;
let currentToken: string;

export async function activate(ctx: vscode.ExtensionContext): Promise<void> {
  logger = vscode.window.createOutputChannel("STAR Launchpad Bridge");
  ctx.subscriptions.push(logger);

  currentToken = await getOrCreateToken(ctx);

  server = new BridgeServer(() => currentToken, logger);
  ctx.subscriptions.push(server);

  statusBar = vscode.window.createStatusBarItem(
    vscode.StatusBarAlignment.Right,
    100
  );
  statusBar.command = "launchpadBridge.showStatus";
  statusBar.text = "$(plug) Bridge";
  statusBar.tooltip = "STAR Launchpad Bridge — starting…";
  statusBar.show();
  ctx.subscriptions.push(statusBar);

  server.onStateChanged((s) => updateStatus(s, server?.clientCount ?? 0));
  server.onClientCountChanged((n) => {
    if (server) {
      const port = server.currentPort;
      updateStatus(
        port
          ? { kind: "listening", port, clients: n }
          : { kind: "stopped" },
        n
      );
    }
  });

  const range = readPortRange();
  try {
    await server.start(range);
  } catch (err) {
    const msg = err instanceof Error ? err.message : String(err);
    logger.appendLine(`Failed to start bridge: ${msg}`);
    vscode.window.showErrorMessage(`STAR Launchpad Bridge: ${msg}`);
  }

  ctx.subscriptions.push(
    vscode.commands.registerCommand("launchpadBridge.pair", () =>
      openPairingPage(ctx)
    ),
    vscode.commands.registerCommand("launchpadBridge.resetToken", async () => {
      currentToken = await rotateToken(ctx);
      server?.invalidateAllSessions("token rotated");
      vscode.window.showInformationMessage(
        "STAR Launchpad Bridge: auth token rotated; existing sessions closed."
      );
    }),
    vscode.commands.registerCommand("launchpadBridge.showStatus", () =>
      showStatusMessage()
    ),
    vscode.commands.registerCommand("launchpadBridge.copyToken", async () => {
      await vscode.env.clipboard.writeText(currentToken);
      vscode.window.showInformationMessage(
        "STAR Launchpad Bridge: auth token copied to clipboard."
      );
    })
  );

  ctx.subscriptions.push(
    vscode.workspace.onDidSaveTextDocument((doc) => {
      server?.broadcast({
        type: "editor_event",
        kind: "file_saved",
        payload: {
          uri: doc.uri.toString(),
          languageId: doc.languageId,
          lineCount: doc.lineCount,
        },
      });
    })
  );

  // First-run pairing nudge (best-effort; users can ignore).
  const firstRun = ctx.globalState.get<boolean>("launchpadBridge.firstRunDone") !== true;
  if (firstRun) {
    await ctx.globalState.update("launchpadBridge.firstRunDone", true);
    void promptFirstRunPairing(ctx);
  }
}

export function deactivate(): void {
  server?.dispose();
  server = null;
  statusBar?.dispose();
  statusBar = null;
}

function readPortRange(): [number, number] {
  const cfg = vscode.workspace.getConfiguration("launchpadBridge");
  const raw = cfg.get<number[]>("portRange", [7777, 7790]);
  if (
    Array.isArray(raw) &&
    raw.length === 2 &&
    Number.isInteger(raw[0]) &&
    Number.isInteger(raw[1])
  ) {
    return [raw[0], raw[1]];
  }
  return [7777, 7790];
}

function updateStatus(state: BridgeState, clients: number): void {
  if (!statusBar) return;
  switch (state.kind) {
    case "listening":
      statusBar.text = clients > 0 ? "$(check) Bridge" : "$(plug) Bridge";
      statusBar.tooltip =
        `STAR Launchpad Bridge — listening on 127.0.0.1:${state.port}\n` +
        `${clients} authed client${clients === 1 ? "" : "s"}.\n` +
        "Click for status.";
      break;
    case "error":
      statusBar.text = "$(warning) Bridge";
      statusBar.tooltip = `STAR Launchpad Bridge error: ${state.message}`;
      break;
    case "stopped":
    default:
      statusBar.text = "$(circle-slash) Bridge";
      statusBar.tooltip = "STAR Launchpad Bridge is stopped.";
      break;
  }
}

function showStatusMessage(): void {
  if (!server) {
    vscode.window.showInformationMessage("STAR Launchpad Bridge: not initialized.");
    return;
  }
  const port = server.currentPort;
  const clients = server.clientCount;
  if (!port) {
    vscode.window.showWarningMessage(
      "STAR Launchpad Bridge: not listening (check the output channel)."
    );
    return;
  }
  void vscode.window.showInformationMessage(
    `STAR Launchpad Bridge listening on 127.0.0.1:${port}. ${clients} authed client(s).`,
    "Open pairing page",
    "Copy token"
  ).then((pick) => {
    if (pick === "Open pairing page") {
      void vscode.commands.executeCommand("launchpadBridge.pair");
    } else if (pick === "Copy token") {
      void vscode.commands.executeCommand("launchpadBridge.copyToken");
    }
  });
}

async function openPairingPage(ctx: vscode.ExtensionContext): Promise<void> {
  const port = server?.currentPort;
  if (!port) {
    vscode.window.showErrorMessage(
      "STAR Launchpad Bridge is not listening; cannot pair."
    );
    return;
  }
  const base = vscode.workspace
    .getConfiguration("launchpadBridge")
    .get<string>("spaUrl", "http://127.0.0.1:8756/launchpad/");
  const url = buildPairUrl(base, port, await getOrCreateToken(ctx));
  await vscode.env.openExternal(vscode.Uri.parse(url));
}

function buildPairUrl(base: string, port: number, token: string): string {
  const trimmed = base.replace(/#.*$/, "");
  const sep = trimmed.includes("?") ? "&" : "?";
  return (
    `${trimmed}${sep}pair=1&port=${port}&token=${encodeURIComponent(token)}` +
    `&label=${encodeURIComponent(vscode.env.appName ?? "editor")}`
  );
}

async function promptFirstRunPairing(ctx: vscode.ExtensionContext): Promise<void> {
  const choice = await vscode.window.showInformationMessage(
    "STAR Launchpad Bridge is running. Pair with the Launchpad now?",
    "Pair now",
    "Copy token",
    "Later"
  );
  if (choice === "Pair now") {
    await openPairingPage(ctx);
  } else if (choice === "Copy token") {
    await vscode.commands.executeCommand("launchpadBridge.copyToken");
  }
}
