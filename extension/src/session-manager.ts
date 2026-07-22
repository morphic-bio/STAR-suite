import * as vscode from "vscode";
import type { WebSocket } from "ws";
import { buildPrompt } from "./prompt";
import type {
  ChatMessage,
  FormData,
  ServerMessage,
  SubmitMode,
} from "./protocol";

/**
 * Accept any registered chat-model provider. We previously constrained this to
 * `{ vendor: "copilot" }` but that locked out Continue, Cline, and
 * provider-neutral LM extensions that register via `vscode.lm`.
 */
const MODEL_PREF: vscode.LanguageModelChatSelector = {};

/** Per-request state kept so we can cancel in-flight LM calls by id. */
interface Session {
  source: vscode.CancellationTokenSource;
  socket: WebSocket;
  mode: SubmitMode;
  history: ChatMessage[];
  assistantAccum: string;
}

/**
 * Owns the mapping from request `id` → active LM request. Extracted from the
 * server so it can be unit-tested with a mocked `vscode.lm` if needed.
 */
export class SessionManager {
  private readonly sessions = new Map<string, Session>();

  constructor(private readonly log: vscode.OutputChannel) {}

  /** Start a new LM request for `id`; the response streams back via `send`. */
  async submit(args: {
    id: string;
    mode: SubmitMode;
    formData: FormData;
    context: ChatMessage[];
    socket: WebSocket;
    send: (msg: ServerMessage) => void;
  }): Promise<void> {
    const { id, mode, formData, context, socket, send } = args;
    if (this.sessions.has(id)) {
      send({ type: "error", id, message: `Duplicate request id: ${id}` });
      return;
    }

    const source = new vscode.CancellationTokenSource();
    const history = [...(context || [])];
    const session: Session = { source, socket, mode, history, assistantAccum: "" };
    this.sessions.set(id, session);
    this.log.appendLine(`Session ${id} started (${mode}); ${history.length} ctx msgs.`);

    send({ type: "ack", id });

    const { user, system } = buildPrompt(mode, formData, history);

    try {
      await this.streamResponse({
        id,
        session,
        userPrompt: user,
        systemPrompt: system,
        send,
      });
    } finally {
      this.sessions.delete(id);
      source.dispose();
    }
  }

  /** Continue an existing thread with a free-text user message. */
  async followup(args: {
    id: string;
    message: string;
    socket: WebSocket;
    send: (msg: ServerMessage) => void;
  }): Promise<void> {
    const { id, message, socket, send } = args;
    if (this.sessions.has(id)) {
      send({
        type: "error",
        id,
        message: "Previous request is still streaming; cancel it first.",
      });
      return;
    }
    const source = new vscode.CancellationTokenSource();
    const session: Session = {
      source,
      socket,
      mode: "refine",
      history: [{ role: "user", content: message }],
      assistantAccum: "",
    };
    this.sessions.set(id, session);

    send({ type: "ack", id });

    try {
      await this.streamResponse({
        id,
        session,
        userPrompt: message,
        systemPrompt:
          "You are continuing a STAR Launchpad bridge conversation. Answer briefly " +
          "and stay focused on the user's follow-up question.",
        send,
      });
    } finally {
      this.sessions.delete(id);
      source.dispose();
    }
  }

  /** Cancel the in-flight request matching `id`, if any. */
  cancel(id: string): boolean {
    const s = this.sessions.get(id);
    if (!s) return false;
    s.source.cancel();
    this.log.appendLine(`Session ${id} cancelled.`);
    return true;
  }

  /** Cancel every in-flight request (e.g. on extension deactivation). */
  cancelAll(): void {
    for (const s of this.sessions.values()) {
      s.source.cancel();
    }
    this.sessions.clear();
  }

  // --- internals ---------------------------------------------------------

  private async streamResponse(args: {
    id: string;
    session: Session;
    userPrompt: string;
    systemPrompt: string;
    send: (msg: ServerMessage) => void;
  }): Promise<void> {
    const { id, session, userPrompt, systemPrompt, send } = args;
    const token = session.source.token;

    let model: vscode.LanguageModelChat | undefined;
    try {
      const models = await vscode.lm.selectChatModels(MODEL_PREF);
      model = models[0];
    } catch (err) {
      send({
        type: "error",
        id,
        message: `Could not pick an LM model: ${describeError(err)}`,
      });
      return;
    }
    if (!model) {
      send({
        type: "error",
        id,
        message:
          "No chat model is available in this editor. Install and sign into a " +
          "provider that registers with vscode.lm (e.g. GitHub Copilot, Continue).",
      });
      return;
    }

    const messages: vscode.LanguageModelChatMessage[] = [
      vscode.LanguageModelChatMessage.User(systemPrompt),
      vscode.LanguageModelChatMessage.User(userPrompt),
    ];

    let response: vscode.LanguageModelChatResponse;
    try {
      response = await model.sendRequest(messages, {}, token);
    } catch (err) {
      if (token.isCancellationRequested) {
        send({ type: "error", id, message: "cancelled" });
      } else {
        send({
          type: "error",
          id,
          message: `LM request failed: ${describeError(err)}`,
        });
      }
      return;
    }

    try {
      for await (const fragment of response.text) {
        if (token.isCancellationRequested) {
          send({ type: "error", id, message: "cancelled" });
          return;
        }
        if (typeof fragment === "string" && fragment.length > 0) {
          session.assistantAccum += fragment;
          send({ type: "chunk", id, text: fragment });
        }
      }
      send({ type: "complete", id });
    } catch (err) {
      if (token.isCancellationRequested) {
        send({ type: "error", id, message: "cancelled" });
        return;
      }
      send({
        type: "error",
        id,
        message: `LM stream error: ${describeError(err)}`,
      });
    }
  }
}

function describeError(err: unknown): string {
  if (err instanceof Error) return err.message;
  try {
    return JSON.stringify(err);
  } catch {
    return String(err);
  }
}
