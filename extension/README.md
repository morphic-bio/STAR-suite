# STAR Launchpad Bridge

VS Code / Cursor extension that pairs with the STAR Launchpad SPA over a
loopback WebSocket so the Launchpad can submit prompts to the editor's chat
and receive streamed responses. This is **phase 1** of the rollout — pairing
and auth only. LM (Copilot Chat / `vscode.lm`) integration ships in phase 2.

## Status

| Feature | Phase 1 | Phase 2 | Phase 3 |
| --- | --- | --- | --- |
| Token auth, loopback WebSocket, pairing URL | ✅ | | |
| `submit` stub that round-trips `ack` → `chunk` → `complete` | ✅ | | |
| `vscode.lm` streaming for refine / research | | ⏳ | |
| `cancel`, `followup` | stub only | ⏳ | |
| Status bar, multi-editor labelling, `editor_event` fanout | ✅ (basic) | | polish |

## Install (developer build)

Requires Node ≥ 18 and VS Code ≥ 1.85.

```bash
cd extension
npm install
npm run compile       # writes out/*.js
```

To package an installable `.vsix`:

```bash
npm run package       # uses vsce (install globally if needed: npm i -g @vscode/vsce)
code --install-extension star-launchpad-bridge-*.vsix
```

To iterate without packaging, open the `extension/` folder in VS Code and
press **F5** — VS Code launches a second "Extension Development Host"
window with the bridge active.

## Settings (`settings.json`)

| Key | Default | Notes |
| --- | --- | --- |
| `launchpadBridge.spaUrl` | `http://127.0.0.1:8756/launchpad/` | Pairing query params are appended to this URL. |
| `launchpadBridge.portRange` | `[7777, 7790]` | Inclusive `[low, high]` TCP range on 127.0.0.1 the bridge may bind to. |

## Commands

- **STAR Launchpad Bridge: Open pairing page** — opens the SPA with
  `?pair=1&port=…&token=…&label=…` so the browser tab can save the editor.
- **STAR Launchpad Bridge: Reset auth token** — rotates the token in
  `globalState` and closes every open session.
- **STAR Launchpad Bridge: Copy auth token** — copies the token to the
  clipboard for manual pairing.
- **STAR Launchpad Bridge: Show status** — info message with port,
  client count, and quick actions.

## Security notes

- The WebSocket server binds to `127.0.0.1` only — it is **not** reachable
  from other machines.
- Clients must send `{ "type": "auth", "token": "…" }` within 2 seconds or
  the server closes the socket.
- The token is 16 random bytes (32 hex chars), stored in VS Code's
  `globalState`, and never logged. Rotate it with
  `launchpadBridge.resetToken` if you suspect it was shared.
- Any process on the same host with the token can connect. This is fine for
  a single-user dev setup; do **not** treat it as a trust boundary against
  other users on a shared machine.

## Protocol (phase 1)

All messages are JSON. See `src/protocol.ts` for canonical types.

```
Client → Server: auth, submit, followup, cancel
Server → Client: auth_ok, auth_fail, ack, chunk, complete, error, editor_event
```

Every `submit` / `followup` carries a UUID `id`; the server echoes that id
on every `ack`, `chunk`, `complete`, and `error` so the SPA can route
streams to the right output panel.

## What phase 2 adds

- Replace the stubbed `submit` handler with `vscode.lm.selectChatModels` +
  `sendRequest`, streaming fragments through `chunk` messages.
- Implement `cancel` via a per-id `CancellationTokenSource` map.
- Add an Anthropic-API fallback behind a setting (optional; skipped for v1
  unless requested).
