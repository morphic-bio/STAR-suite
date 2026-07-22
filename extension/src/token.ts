import * as crypto from "crypto";
import * as vscode from "vscode";

const TOKEN_KEY = "launchpadBridge.token";

/** Ensure a token exists in globalState; return it. */
export async function getOrCreateToken(ctx: vscode.ExtensionContext): Promise<string> {
  const existing = ctx.globalState.get<string>(TOKEN_KEY);
  if (existing && /^[0-9a-f]{32}$/.test(existing)) {
    return existing;
  }
  const fresh = crypto.randomBytes(16).toString("hex");
  await ctx.globalState.update(TOKEN_KEY, fresh);
  return fresh;
}

/** Replace the stored token with a fresh one. */
export async function rotateToken(ctx: vscode.ExtensionContext): Promise<string> {
  const fresh = crypto.randomBytes(16).toString("hex");
  await ctx.globalState.update(TOKEN_KEY, fresh);
  return fresh;
}

/** Constant-time compare to avoid leaking timing info on mismatch. */
export function tokenMatches(expected: string, candidate: unknown): boolean {
  if (typeof candidate !== "string") return false;
  if (candidate.length !== expected.length) return false;
  try {
    return crypto.timingSafeEqual(
      Buffer.from(expected, "utf8"),
      Buffer.from(candidate, "utf8")
    );
  } catch {
    return false;
  }
}
