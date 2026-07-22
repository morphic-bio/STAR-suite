/**
 * Bridge message protocol shared by the VS Code extension (WebSocket server)
 * and the STAR Launchpad SPA (WebSocket client). Keep this file in sync with
 * the JSON shapes used by mcp_server/launchpad/static/bridge.js.
 */

/** Chat context item sent alongside a submit request. */
export interface ChatMessage {
  role: "user" | "assistant" | "system";
  content: string;
}

/** Generic form payload; the Launchpad defines its own internal shape. */
export type FormData = Record<string, unknown>;

export type SubmitMode = "refine" | "research";

// SPA -> extension ----------------------------------------------------------
export interface AuthMsg {
  type: "auth";
  token: string;
}

export interface SubmitMsg {
  type: "submit";
  id: string;
  mode: SubmitMode;
  formData: FormData;
  context: ChatMessage[];
}

export interface FollowupMsg {
  type: "followup";
  id: string;
  message: string;
}

export interface CancelMsg {
  type: "cancel";
  id: string;
}

export type ClientMessage = AuthMsg | SubmitMsg | FollowupMsg | CancelMsg;

// extension -> SPA ----------------------------------------------------------
export interface AuthOkMsg {
  type: "auth_ok";
  /** VS Code / Cursor app name (e.g. "Visual Studio Code" / "Cursor"). */
  app: string;
  /** Workspace folder name, or empty string. */
  workspace: string;
  /** os.hostname() — used by the SPA to label the editor. */
  host: string;
  /** Protocol version so the SPA can gate new features. */
  protocol: number;
}

export interface AuthFailMsg {
  type: "auth_fail";
  /** Human-readable reason (e.g. "token mismatch", "auth timeout"). */
  reason: string;
}

export interface AckMsg {
  type: "ack";
  id: string;
}

export interface ChunkMsg {
  type: "chunk";
  id: string;
  text: string;
}

export interface CompleteMsg {
  type: "complete";
  id: string;
}

export interface ErrorMsg {
  type: "error";
  id: string;
  message: string;
}

export interface EditorEventMsg {
  type: "editor_event";
  kind: "file_saved" | "selection_changed";
  payload: Record<string, unknown>;
}

export type ServerMessage =
  | AuthOkMsg
  | AuthFailMsg
  | AckMsg
  | ChunkMsg
  | CompleteMsg
  | ErrorMsg
  | EditorEventMsg;

/** Bump this when incompatible message-shape changes land. */
export const PROTOCOL_VERSION = 1;
