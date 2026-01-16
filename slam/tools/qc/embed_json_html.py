#!/usr/bin/env python3
"""
Make STAR-Flex QC HTML self-contained by embedding the adjacent JSON and removing fetch().

This is useful for opening reports directly in Firefox via file://, where fetch() is blocked.

Usage:
  tools/qc/embed_json_html.py --html in.html --json in.json --out out.html

Notes:
- This script does NOT bundle Plotly JS. The HTML still loads Plotly from the CDN.
- Newer STAR-Flex builds embed JSON by default; this script is for older outputs.
"""

import argparse
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--html", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    html_path = Path(args.html)
    json_path = Path(args.json)
    out_path = Path(args.out)

    html = html_path.read_text()
    json_text = json_path.read_text()

    if 'type="application/json" id="qc_json"' in html:
        # Already embedded; just copy through.
        out_path.write_text(html)
        return 0

    # Insert embedded JSON before the first <script> that begins with fetch(
    fetch_idx = html.find("fetch(")
    if fetch_idx == -1:
        raise SystemExit("Could not find fetch() in HTML; nothing to patch")
    script_idx = html.rfind("<script", 0, fetch_idx)
    if script_idx == -1:
        raise SystemExit("Could not find <script> before fetch()")

    embed = (
        '  <script type="application/json" id="qc_json">\n'
        + json_text
        + ("" if json_text.endswith("\n") else "\n")
        + "  </script>\n"
    )

    html2 = html[:script_idx] + embed + html[script_idx:]

    # Now remove fetch wrapper in the first script that contains it.
    start = html2.find("fetch(")
    if start == -1:
        out_path.write_text(html2)
        return 0

    # Replace leading fetch chain with JSON.parse.
    then_marker = "then(data => {"
    then_idx = html2.find(then_marker, start)
    if then_idx == -1:
        then_marker = "then(data=>{"
        then_idx = html2.find(then_marker, start)
    if then_idx == -1:
        raise SystemExit("Could not find then(data => {) in fetch chain")

    # Position after the opening brace of then(data => { ... })
    body_start = then_idx + len(then_marker)

    prefix = (
        "try {\n"
        "  if (typeof Plotly === 'undefined') { throw new Error('Plotly failed to load. If you are offline or the CDN is blocked, the plots cannot render.'); }\n"
        "  const data = JSON.parse(document.getElementById('qc_json').textContent);\n"
    )

    html2 = html2[:start] + prefix + html2[body_start:]

    # Remove the terminal `});` that closes the fetch chain (just before </script>)
    script_close = html2.find("</script>", start)
    if script_close == -1:
        raise SystemExit("Could not find </script> after fetch()")
    tail = html2.rfind("});", start, script_close)
    if tail == -1:
        raise SystemExit("Could not find closing }); for fetch chain")

    suffix = (
        "\n} catch (e) {\n"
        "  const el = document.getElementById('summary');\n"
        "  if (el) el.innerHTML = '<div style=\"color:#b00020;font-weight:bold\">Failed to render plots: ' + (e && e.message ? e.message : String(e)) + '</div>';\n"
        "}\n"
    )

    html2 = html2[:tail] + suffix + html2[tail + 3:]

    out_path.write_text(html2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

