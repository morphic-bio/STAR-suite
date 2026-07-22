import type { ChatMessage, FormData, SubmitMode } from "./protocol";

const MAX_CONTEXT_MESSAGES = 10;
const INSTRUCTIONS: Record<SubmitMode, string> = {
  refine:
    "You are reviewing a STAR-suite workflow configuration. For each non-default " +
    "parameter, explain what it does, flag likely issues (e.g. unrealistic values, " +
    "mismatched mate lists), and suggest concrete changes. Keep the reply short and " +
    "structured as a bulleted list grouped by parameter.",
  research:
    "You are assisting a bioinformatician who is about to run a STAR-suite pipeline. " +
    "Using the parameters below, produce (a) a ready-to-run bash snippet that invokes " +
    "STAR with these flags, (b) a short explanation of each stage, and (c) the most " +
    "important caveats for this configuration. Fence the snippet in ```bash.",
};

export interface BuiltPrompt {
  /** Single user-facing prompt string. */
  user: string;
  /** Optional system preamble (vendor-specific; may be prepended). */
  system: string;
}

/**
 * Build a deterministic prompt from form values and the last N chat messages.
 * Kept pure and dep-free so it can be unit-tested and logged.
 */
export function buildPrompt(
  mode: SubmitMode,
  formData: FormData,
  context: ChatMessage[]
): BuiltPrompt {
  const params = formatFormData(formData || {});
  const history = formatContext(context || []);
  const instruction = INSTRUCTIONS[mode] ?? INSTRUCTIONS.refine;

  const user = [
    `# Task`,
    mode,
    "",
    `# User inputs`,
    params,
    "",
    `# Prior conversation`,
    history,
    "",
    `# Instruction`,
    instruction,
  ].join("\n");

  const system =
    "You are the AI assistant paired with STAR Launchpad, a UI for running " +
    "STAR alignment workflows. Always answer in the context of the STAR CLI " +
    "and the provided parameters. Never invent parameter values the user did " +
    "not supply.";

  return { user, system };
}

function formatFormData(form: FormData): string {
  const keys = Object.keys(form);
  if (!keys.length) return "(none)";
  const lines: string[] = [];
  for (const k of keys) {
    const v = form[k];
    lines.push(`- ${k}: ${renderValue(v)}`);
  }
  return lines.join("\n");
}

function formatContext(context: ChatMessage[]): string {
  const trimmed = context.slice(-MAX_CONTEXT_MESSAGES);
  if (!trimmed.length) return "(none)";
  return trimmed
    .map((m) => {
      const role = m.role === "assistant" ? "Assistant" : m.role === "system" ? "System" : "User";
      const text = String(m.content || "").trim();
      return `${role}: ${text}`;
    })
    .join("\n\n");
}

function renderValue(v: unknown): string {
  if (v === null || v === undefined) return "(unset)";
  if (typeof v === "string") return v.length ? v : "(empty)";
  if (typeof v === "number" || typeof v === "boolean") return String(v);
  try {
    return JSON.stringify(v);
  } catch {
    return String(v);
  }
}
