/**
 * STAR Launchpad — schema-driven form; validation/rendering stays server-side.
 * Presentation-only: displayGroups, gated_by / ui_gated_by, constraint help.
 */

/**
 * Resolve JSON API base (e.g. "/launchpad/api" or "/mcp/launchpad/api").
 * Absolute "/launchpad/api/..." breaks when the MCP app is mounted under a path prefix.
 */
function launchpadApiBase() {
  let p = window.location.pathname;
  p = p.replace(/\/+$/, "");
  if (p.endsWith("/index.html")) {
    p = p.slice(0, -"/index.html".length);
  }
  const i = p.lastIndexOf("/launchpad");
  if (i >= 0) {
    return p.slice(0, i + "/launchpad".length) + "/api";
  }
  return "/launchpad/api";
}

/** Core STAR CLI recipes first (order matches typical use: index → align modes). */
const STAR_WORKFLOW_ORDER = [
  "star_genome_generate",
  "star_bulk_pe_batch",
  "star_scrna_solo_droplet",
  "star_flex_fixed_rna",
  "star_perturb_cr_compat",
];

function starWorkflowRank(id) {
  const i = STAR_WORKFLOW_ORDER.indexOf(id);
  return i >= 0 ? i : 999;
}

function sortWorkflowsForLaunchpad(list) {
  return list.slice().sort((a, b) => {
    const aStar = String(a.id || "").startsWith("star_");
    const bStar = String(b.id || "").startsWith("star_");
    if (aStar !== bStar) return aStar ? -1 : 1;
    if (aStar && bStar) {
      const ra = starWorkflowRank(a.id);
      const rb = starWorkflowRank(b.id);
      if (ra !== rb) return ra - rb;
    }
    return String(a.id || "").localeCompare(String(b.id || ""));
  });
}

function shellQuote(s) {
  return "'" + String(s).replace(/'/g, "'\\''") + "'";
}

/** Minimal annotated Bash used as the Script Lane placeholder (harmless ``cp`` step). */
function scriptLaneDefaultScript() {
  return [
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "#@workflow id=launchpad_script_lane",
    "#@input genome file required",
    "#@output bam file",
    "#@step id=copy tool=cp",
    "#@uses genome",
    "#@produces bam",
    'cp "${genome}" "${bam}"',
    "",
  ].join("\n");
}

function launchpadApp() {
  return {
    activeTab: "workflows",
    /** Full list from API (sorted); `workflows` is filtered for the dropdown. */
    allWorkflows: [],
    workflows: [],
    /** When false, recipe list is only `star_*` CLI recipes. */
    includeTestWorkflows: false,
    workflowId: "",
    schema: null,
    describe: null,
    displayGroups: [],
    params: {},
    focusParam: null,
    loading: false,
    validateResult: null,
    renderResult: null,
    launchResult: null,
    launchSupported: false,
    serverPathCheckSupported: false,
    scriptLaneUtilsReady: false,
    scriptLaneViabilityInputsSupported: false,
    scriptLaneExecuteSupported: false,
    /** True when this browser is on the same host as the server (loopback). */
    trustedLocal: false,
    /** Whether the server exposes /browse (true when trusted_roots is non-empty). */
    browseSupported: false,
    uploadSupported: false,
    /** Resolved trusted_roots, shown at the top of the server file picker. */
    browseRoots: [],
    checkPaths: false,
    lastValidationCheckPaths: null,
    httpError: "",
    notice: "",
    showAbout: false,

    // --- Server file picker modal state ----------------------------------
    showFilePicker: false,
    /** Name of the param we'll write the chosen path into. */
    fpTargetParam: null,
    /** 'file' | 'directory' | 'string_list' — decides what is selectable. */
    fpMode: "file",
    /** Current directory being viewed; empty string = virtual roots list. */
    fpPath: "",
    fpParent: null,
    fpEntries: [],
    /** {name, kind, size, mtime} of entries the user has clicked. */
    fpSelected: [],
    fpTruncated: false,
    fpLoading: false,
    fpError: "",
    /** File input used by the local-file upload flow. */
    uploadBusy: false,
    uploadTargetParam: null,

    slScript: scriptLaneDefaultScript(),
    slInputJson: "{}",
    slWorkdir: "",
    slTranslateResult: null,
    slViabilityResult: null,
    slExecuteResult: null,
    slIr: null,
    slNormalizedScript: null,

    /** Drop IR and downstream Script Lane state (script edited, translate failed, etc.). */
    slInvalidateScriptLaneDerivedState() {
      this.slIr = null;
      this.slNormalizedScript = null;
      this.slTranslateResult = null;
      this.slViabilityResult = null;
      this.slExecuteResult = null;
    },

    selectedParam() {
      if (!this.focusParam) return null;
      return this.paramByName(this.focusParam);
    },

    paramByName(name) {
      if (!this.schema?.parameters) return null;
      return this.schema.parameters.find((p) => p.name === name) || null;
    },

    /** Native tooltip text (schema description + optional env / choices). */
    paramTooltip(pname) {
      const p = this.paramByName(pname);
      if (!p) return "";
      const lines = [];
      const desc = String(p.description || "").trim();
      if (desc) lines.push(desc);
      if (p.env_var) lines.push(`Environment variable: ${p.env_var}`);
      if (p.choices && p.choices.length)
        lines.push(`Allowed values: ${p.choices.join(", ")}`);
      if (lines.length) return lines.join("\n\n");
      if (p.cli_flag) return `Flag ${p.cli_flag} (${p.name}).`;
      return pname;
    },

    formatMetaValue(value) {
      if (value === null || value === undefined) return "";
      if (typeof value === "boolean") return value ? "true" : "false";
      if (Array.isArray(value)) return value.join(", ");
      if (typeof value === "object") return JSON.stringify(value);
      return String(value);
    },

    paramRangeText(p) {
      if (!p) return "";
      const lo = p.min_value;
      const hi = p.max_value;
      if (lo === null || lo === undefined) {
        if (hi === null || hi === undefined) return "";
        return `<= ${this.formatMetaValue(hi)}`;
      }
      if (hi === null || hi === undefined) return `>= ${this.formatMetaValue(lo)}`;
      return `${this.formatMetaValue(lo)} to ${this.formatMetaValue(hi)}`;
    },

    paramAliasesText(p) {
      if (!p?.aliases?.length) return "";
      return p.aliases.join(", ");
    },

    normalizedParamsText() {
      const obj = this.validateResult?.normalized_params;
      if (!obj) return "";
      return JSON.stringify(obj, null, 2);
    },

    validationModeMessage() {
      if (this.lastValidationCheckPaths === true) {
        return "Validation included file and directory existence checks on the server host.";
      }
      if (this.lastValidationCheckPaths === false) {
        return "Validation checked schema and types only. File and directory existence were not checked.";
      }
      return "";
    },

    paramHasRichTip(pname) {
      const p = this.paramByName(pname);
      if (!p) return false;
      if (String(p.description || "").trim()) return true;
      if (String(p.help || "").trim()) return true;
      if (p.env_var) return true;
      if (p.choices && p.choices.length) return true;
      return false;
    },

    /** True when a boolean gate parameter is on (e.g. star_only). */
    gateIsActive(gateName) {
      if (!gateName) return false;
      const v = this.params[gateName];
      if (typeof v === "boolean") return v;
      if (v === "" || v == null) return false;
      if (typeof v === "number") return v !== 0;
      const s = String(v).toLowerCase();
      return s === "true" || s === "1" || s === "yes";
    },

    groupHidden(group) {
      return !!(group && group.gated_by && this.gateIsActive(group.gated_by));
    },

    paramRowVisible(pname) {
      const p = this.paramByName(pname);
      if (!p) return false;
      if (p.ui_gated_by && this.gateIsActive(p.ui_gated_by)) return false;
      return true;
    },

    /** True if the param would appear in the form (visible group + no ui gate). */
    paramVisibleInForm(pname) {
      if (!this.paramRowVisible(pname)) return false;
      for (const g of this.displayGroups || []) {
        if (!(g.parameters || []).includes(pname)) continue;
        if (this.groupHidden(g)) return false;
      }
      return true;
    },

    firstVisibleParamName() {
      for (const g of this.displayGroups || []) {
        if (this.groupHidden(g)) continue;
        for (const pname of g.parameters || []) {
          if (this.paramRowVisible(pname) && this.paramByName(pname)) return pname;
        }
      }
      return null;
    },

    /** Clear or move focus when gates hide the current field (help pane stays in sync). */
    ensureFocusVisible() {
      if (!this.focusParam) return;
      if (this.paramVisibleInForm(this.focusParam)) return;
      this.focusParam = this.firstVisibleParamName();
    },

    /** Stages omitted from help when their gate is active (e.g. downstream when star_only). */
    visibleStages() {
      const stages = this.describe?.stages || [];
      return stages.filter((s) => {
        if (!s.gated_by) return true;
        return !this.gateIsActive(s.gated_by);
      });
    },

    /** Constraints that mention the focused parameter. */
    constraintsForFocus() {
      if (
        !this.focusParam ||
        !this.paramVisibleInForm(this.focusParam) ||
        !this.schema?.constraints
      ) {
        return [];
      }
      return this.schema.constraints.filter(
        (c) => c.params && c.params.includes(this.focusParam)
      );
    },

    rebuildDisplayGroups() {
      if (!this.schema) {
        this.displayGroups = [];
        return;
      }
      const inGroup = new Set();
      for (const g of this.schema.parameter_groups || []) {
        for (const p of g.parameters || []) inGroup.add(p);
      }
      const allNames = (this.schema.parameters || []).map((p) => p.name);
      const other = allNames.filter((n) => !inGroup.has(n));
      const base = [...(this.schema.parameter_groups || [])];
      if (other.length) {
        base.push({
          name: "_other",
          title: "Other",
          parameters: other,
          gated_by: null,
        });
      }
      this.displayGroups = base;
      this.ensureFocusVisible();
    },

    renderBlock() {
      if (!this.renderResult) return "";
      const lines = [];
      const env = this.renderResult.env_overrides || {};
      for (const [k, v] of Object.entries(env)) {
        lines.push(`export ${k}=${shellQuote(v)}`);
      }
      if (lines.length) lines.push("");
      lines.push(this.renderResult.shell_preview || "");
      return lines.join("\n");
    },

    async copyCmd() {
      const t = this.renderBlock();
      try {
        await navigator.clipboard.writeText(t);
      } catch (e) {
        this.httpError = "Clipboard unavailable: " + (e && e.message ? e.message : String(e));
      }
    },

    defaultParams() {
      const out = {};
      for (const p of this.schema.parameters) {
        const d = p.default;
        if (d === null || d === undefined) {
          if (p.type === "bool") out[p.name] = false;
          else out[p.name] = "";
        } else {
          out[p.name] = d;
        }
      }
      if (this.workflowId === "ucsf_star_suite_production") {
        out.star_only = true;
      }
      return out;
    },

    collectParams() {
      const out = {};
      for (const p of this.schema.parameters) {
        let v = this.params[p.name];
        if (v === "" || v === null || v === undefined) {
          if (p.type === "bool") {
            out[p.name] = false;
          }
          continue;
        }
        if (p.type === "bool") {
          out[p.name] = !!v;
        } else if (p.type === "int") {
          const n = parseInt(String(v), 10);
          if (!Number.isNaN(n)) out[p.name] = n;
        } else if (p.type === "float") {
          const n = parseFloat(String(v));
          if (!Number.isNaN(n)) out[p.name] = n;
        } else {
          out[p.name] = v;
        }
      }
      return out;
    },

    /** Merge saved values onto schema defaults (drops unknown keys). */
    applySavedParamValues(raw) {
      const out = this.defaultParams();
      if (!raw || typeof raw !== "object") return out;
      for (const p of this.schema.parameters) {
        if (!Object.prototype.hasOwnProperty.call(raw, p.name)) continue;
        let v = raw[p.name];
        if (p.type === "bool") {
          out[p.name] = !!v;
        } else if (p.type === "int") {
          const n = parseInt(String(v), 10);
          out[p.name] = Number.isNaN(n) ? out[p.name] : n;
        } else if (p.type === "float") {
          const n = parseFloat(String(v));
          out[p.name] = Number.isNaN(n) ? out[p.name] : n;
        } else if (v === null || v === undefined) {
          continue;
        } else {
          out[p.name] = v;
        }
      }
      return out;
    },

    /** Pick a local JSON file and merge matching parameter keys into the form (client-only). */
    loadParametersFromFile(event) {
      this.httpError = "";
      const el = event.target;
      const file = el.files && el.files[0];
      el.value = "";
      if (!file || !this.schema?.parameters) return;
      const reader = new FileReader();
      reader.onload = () => {
        try {
          const data = JSON.parse(String(reader.result || ""));
          const rawParams =
            data &&
            typeof data === "object" &&
            !Array.isArray(data) &&
            data.params &&
            typeof data.params === "object" &&
            !Array.isArray(data.params)
              ? data.params
              : data;
          if (!rawParams || typeof rawParams !== "object" || Array.isArray(rawParams)) {
            this.httpError = "JSON must be an object with parameter keys or a { params: { … } } wrapper.";
            return;
          }
          const names = new Set(this.schema.parameters.map((p) => p.name));
          let hit = 0;
          for (const k of Object.keys(rawParams)) {
            if (names.has(k)) hit++;
          }
          if (hit === 0) {
            this.httpError = "No keys in file match this recipe’s parameters.";
            return;
          }
          if (data.workflow_id && data.workflow_id !== this.workflowId) {
            this.httpError =
              `File is for recipe "${data.workflow_id}"; applied overlapping fields to "${this.workflowId}".`;
          } else {
            this.httpError = "";
          }
          this.params = this.applySavedParamValues(rawParams);
          this.validateResult = null;
          this.renderResult = null;
          this.launchResult = null;
          this.rebuildDisplayGroups();
          this.focusParam = this.firstVisibleParamName();
        } catch (e) {
          this.httpError = "Invalid JSON: " + (e && e.message ? e.message : String(e));
        }
      };
      reader.onerror = () => {
        this.httpError = "Could not read file.";
      };
      reader.readAsText(file);
    },

    /** Download current form values as JSON (client-only). */
    saveParametersToFile() {
      this.httpError = "";
      if (!this.schema?.parameters) return;
      const snapshot = {};
      for (const p of this.schema.parameters) {
        snapshot[p.name] = this.params[p.name];
      }
      const payload = {
        workflow_id: this.workflowId,
        exported_at: new Date().toISOString(),
        params: snapshot,
      };
      const blob = new Blob([JSON.stringify(payload, null, 2)], {
        type: "application/json",
      });
      const a = document.createElement("a");
      const url = URL.createObjectURL(blob);
      a.href = url;
      a.download = `launchpad-${this.workflowId}-params.json`;
      a.click();
      URL.revokeObjectURL(url);
    },

    // --- File picker / upload ------------------------------------------
    /** True if a parameter should get a server-browse button. */
    paramSupportsServerBrowse(p) {
      if (!p) return false;
      return ["file", "directory", "string_list"].includes(p.type);
    },

    /** True if a parameter should get a local-file upload button. */
    paramSupportsLocalUpload(p) {
      if (!p || this.trustedLocal) return false;
      if (!this.uploadSupported) return false;
      // Uploads produce a single server path; only meaningful for file / list-of-files.
      return p.type === "file" || p.type === "string_list";
    },

    /** True if both explorers would appear for this param (for layout). */
    paramHasBothExplorers(p) {
      return this.paramSupportsServerBrowse(p) && this.paramSupportsLocalUpload(p);
    },

    formatFileSize(bytes) {
      if (bytes == null) return "";
      const n = Number(bytes);
      if (!Number.isFinite(n) || n < 0) return "";
      if (n < 1024) return `${n} B`;
      const units = ["KB", "MB", "GB", "TB"];
      let v = n / 1024;
      let i = 0;
      while (v >= 1024 && i < units.length - 1) {
        v /= 1024;
        i++;
      }
      const rounded = v >= 10 ? Math.round(v) : Math.round(v * 10) / 10;
      return `${rounded} ${units[i]}`;
    },

    /** Open the server file picker modal for this param. */
    openFilePicker(pname) {
      const p = this.paramByName(pname);
      if (!p) return;
      if (!this.browseSupported) {
        this.httpError = "Server filesystem browsing is not available on this server.";
        return;
      }
      this.fpTargetParam = pname;
      this.fpMode = p.type;
      this.fpSelected = [];
      this.fpError = "";
      this.fpTruncated = false;
      this.fpEntries = [];
      this.fpPath = "";
      this.fpParent = null;
      this.showFilePicker = true;
      // Seed with a sensible starting directory: current value if it exists, else roots.
      const current = String(this.params[pname] || "").trim();
      // For string_list pick the first entry as starting dir.
      const seed = current.split(",")[0].trim();
      const startDir = seed ? this._dirnameGuess(seed) : "";
      void this.fpNavigate(startDir);
    },

    closeFilePicker() {
      this.showFilePicker = false;
      this.fpTargetParam = null;
      this.fpSelected = [];
      this.fpEntries = [];
      this.fpError = "";
    },

    /** Heuristic: strip the last path component for seeding browse-at-current-value. */
    _dirnameGuess(p) {
      if (!p) return "";
      const s = String(p).replace(/\/+$/, "");
      const i = s.lastIndexOf("/");
      if (i < 0) return "";
      return i === 0 ? "/" : s.slice(0, i);
    },

    async fpNavigate(path) {
      this.fpLoading = true;
      this.fpError = "";
      try {
        const api = launchpadApiBase();
        const url = path
          ? `${api}/browse?path=${encodeURIComponent(path)}`
          : `${api}/browse`;
        const r = await fetch(url);
        const data = await r.json().catch(() => ({}));
        if (!r.ok) {
          // Fall back to root list on forbidden/missing so picker stays usable.
          this.fpError = data.message || r.statusText;
          if (path) await this.fpNavigate("");
          return;
        }
        this.fpPath = data.is_root_list ? "" : (data.path || "");
        this.fpParent = data.parent || null;
        this.fpEntries = Array.isArray(data.entries) ? data.entries : [];
        this.fpTruncated = !!data.truncated;
      } catch (e) {
        this.fpError = e && e.message ? e.message : String(e);
      } finally {
        this.fpLoading = false;
      }
    },

    /** Double-click / Enter on an entry: navigate into dirs, toggle selection on files. */
    fpActivate(entry) {
      if (!entry) return;
      if (entry.kind === "dir") {
        const target = this.fpPath === "" ? entry.name : this._joinPath(this.fpPath, entry.name);
        void this.fpNavigate(target);
        return;
      }
      this.fpToggleSelect(entry);
    },

    /** Single click: select; semantics depend on param type. */
    fpToggleSelect(entry) {
      if (!entry) return;
      const selectable = this._entrySelectable(entry);
      if (!selectable) return;
      const fullPath = this._entryFullPath(entry);
      const multi = this.fpMode === "string_list";
      const idx = this.fpSelected.findIndex((s) => s.path === fullPath);
      if (idx >= 0) {
        this.fpSelected.splice(idx, 1);
        return;
      }
      const record = { path: fullPath, kind: entry.kind, name: entry.name };
      if (multi) {
        this.fpSelected.push(record);
      } else {
        this.fpSelected = [record];
      }
    },

    _entrySelectable(entry) {
      if (!entry) return false;
      if (entry.kind === "dir") {
        return this.fpMode === "directory" || this.fpMode === "string_list";
      }
      if (entry.kind === "file") {
        return this.fpMode === "file" || this.fpMode === "string_list";
      }
      return false;
    },

    _entryFullPath(entry) {
      if (!entry) return "";
      if (this.fpPath === "") return entry.name;
      return this._joinPath(this.fpPath, entry.name);
    },

    _joinPath(a, b) {
      if (!a) return b;
      if (a.endsWith("/")) return a + b;
      return a + "/" + b;
    },

    /** "Select this directory" shortcut for directory params. */
    fpPickCurrentDirectory() {
      if (!this.fpPath || this.fpMode !== "directory") return;
      this.fpSelected = [{ path: this.fpPath, kind: "dir", name: this.fpPath }];
      this.fpApply();
    },

    /** Write the selection into the bound parameter and close. */
    fpApply() {
      if (!this.fpTargetParam) return this.closeFilePicker();
      if (!this.fpSelected.length) {
        this.fpError = "Pick at least one entry first.";
        return;
      }
      const paths = this.fpSelected.map((s) => s.path);
      if (this.fpMode === "string_list") {
        this.params[this.fpTargetParam] = paths.join(",");
      } else {
        this.params[this.fpTargetParam] = paths[0];
      }
      this.validateResult = null;
      this.renderResult = null;
      this.launchResult = null;
      this.closeFilePicker();
    },

    fpSelectedPath(entry) {
      return this.fpSelected.some((s) => s.path === this._entryFullPath(entry));
    },

    /** Upload a file chosen on the user's machine; store returned server path. */
    async handleLocalUpload(pname, event) {
      const el = event.target;
      const file = el.files && el.files[0];
      el.value = "";
      if (!file || !pname) return;
      const p = this.paramByName(pname);
      if (!p) return;
      this.uploadTargetParam = pname;
      this.uploadBusy = true;
      this.httpError = "";
      try {
        const fd = new FormData();
        fd.append("file", file, file.name);
        const api = launchpadApiBase();
        const r = await fetch(`${api}/upload`, { method: "POST", body: fd });
        const data = await r.json().catch(() => ({}));
        if (!r.ok || !data.ok) {
          this.httpError = data.message || `Upload failed (${r.status}).`;
          return;
        }
        const serverPath = String(data.path || "");
        if (!serverPath) {
          this.httpError = "Upload succeeded but the server did not return a path.";
          return;
        }
        if (p.type === "string_list") {
          const existing = String(this.params[pname] || "").trim();
          this.params[pname] = existing ? `${existing},${serverPath}` : serverPath;
        } else {
          this.params[pname] = serverPath;
        }
        this.validateResult = null;
        this.renderResult = null;
        this.launchResult = null;
        this.notice = `Uploaded ${data.filename} → ${serverPath}`;
      } catch (e) {
        this.httpError = e && e.message ? e.message : String(e);
      } finally {
        this.uploadBusy = false;
        this.uploadTargetParam = null;
      }
    },

    applyWorkflowFilter() {
      if (this.includeTestWorkflows) {
        this.workflows = this.allWorkflows.slice();
      } else {
        this.workflows = this.allWorkflows.filter((w) =>
          String(w.id || "").startsWith("star_")
        );
      }
    },

    /** After toggling test workflows, keep a valid selection. */
    async onIncludeTestWorkflowsChange() {
      this.applyWorkflowFilter();
      const ids = new Set((this.workflows || []).map((w) => w.id));
      if (ids.has(this.workflowId)) return;
      const star = this.workflows.find((w) => w.id === "star_genome_generate");
      this.workflowId = star
        ? star.id
        : this.workflows.length
          ? this.workflows[0].id
          : "";
      if (this.workflowId) await this.loadWorkflow();
    },

    async init() {
      const api = launchpadApiBase();
      const [wr, cap] = await Promise.all([
        fetch(`${api}/workflows`),
        fetch(`${api}/capabilities`),
      ]);
      if (cap.ok) {
        const c = await cap.json();
        this.launchSupported = !!c.launch_supported;
        this.serverPathCheckSupported = !!c.server_path_check_supported;
        this.scriptLaneUtilsReady = !!c.script_lane_translate_supported;
        this.scriptLaneViabilityInputsSupported = !!c.script_lane_viability_inputs_supported;
        this.scriptLaneExecuteSupported = !!c.script_lane_execute_supported;
        this.trustedLocal = !!c.trusted_local;
        this.browseSupported = !!c.browse_supported;
        this.uploadSupported = !!c.upload_supported;
        this.browseRoots = Array.isArray(c.browse_roots) ? c.browse_roots : [];
        if (!this.serverPathCheckSupported) this.checkPaths = false;
      } else {
        this.launchSupported = false;
        this.serverPathCheckSupported = false;
        this.scriptLaneUtilsReady = false;
        this.scriptLaneViabilityInputsSupported = false;
        this.scriptLaneExecuteSupported = false;
        this.trustedLocal = false;
        this.browseSupported = false;
        this.uploadSupported = false;
        this.browseRoots = [];
        this.checkPaths = false;
      }
      if (!wr.ok) {
        this.httpError = await wr.text();
        return;
      }
      const data = await wr.json();
      const raw = data.workflows || [];
      this.allWorkflows = sortWorkflowsForLaunchpad(raw);
      this.applyWorkflowFilter();
      if (this.workflows.length) {
        const starDefault = this.workflows.find((w) => w.id === "star_genome_generate");
        this.workflowId = starDefault ? starDefault.id : this.workflows[0].id;
        await this.loadWorkflow();
      }
    },

    async loadWorkflow() {
      this.loading = true;
      this.httpError = "";
      this.validateResult = null;
      this.lastValidationCheckPaths = null;
      this.renderResult = null;
      this.launchResult = null;
      try {
        const api = launchpadApiBase();
        const [s, d] = await Promise.all([
          fetch(`${api}/workflows/${encodeURIComponent(this.workflowId)}/schema`),
          fetch(`${api}/workflows/${encodeURIComponent(this.workflowId)}/describe`),
        ]);
        if (!s.ok) {
          this.httpError = (await s.json()).message || (await s.text());
          return;
        }
        if (!d.ok) {
          this.httpError = (await d.json()).message || (await d.text());
          return;
        }
        this.schema = await s.json();
        this.describe = await d.json();
        this.params = this.defaultParams();
        this.rebuildDisplayGroups();
        this.focusParam = this.firstVisibleParamName();
      } finally {
        this.loading = false;
      }
    },

    async doValidate(opts = {}) {
      const checkPaths = !!opts.checkPaths;
      this.loading = true;
      this.httpError = "";
      this.lastValidationCheckPaths = checkPaths;
      try {
        const api = launchpadApiBase();
        const r = await fetch(
          `${api}/workflows/${encodeURIComponent(this.workflowId)}/validate`,
          {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({
              params: this.collectParams(),
              check_paths: checkPaths,
            }),
          }
        );
        const data = await r.json();
        if (!r.ok) {
          this.httpError = data.message || r.statusText;
          return false;
        }
        this.validateResult = data;
        return !!data.valid;
      } catch (e) {
        this.httpError = e && e.message ? e.message : String(e);
        return false;
      } finally {
        this.loading = false;
      }
    },

    async doRender() {
      const valid = await this.doValidate({
        checkPaths: this.serverPathCheckSupported && this.checkPaths,
      });
      if (!valid) return;
      this.loading = true;
      this.httpError = "";
      try {
        const api = launchpadApiBase();
        const r = await fetch(
          `${api}/workflows/${encodeURIComponent(this.workflowId)}/render`,
          {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ params: this.collectParams() }),
          }
        );
        const data = await r.json();
        if (!r.ok) {
          this.httpError = data.message || r.statusText;
          return;
        }
        this.renderResult = data;
        this.launchResult = null;
      } finally {
        this.loading = false;
      }
    },

    async doLaunch() {
      const valid = await this.doValidate({ checkPaths: true });
      if (!valid) return;
      this.loading = true;
      this.httpError = "";
      this.notice = "";
      this.launchResult = null;
      try {
        const api = launchpadApiBase();
        const r = await fetch(
          `${api}/workflows/${encodeURIComponent(this.workflowId)}/launch`,
          {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({
              params: this.collectParams(),
              check_paths: true,
            }),
          }
        );
        const data = await r.json();
        if (r.status === 403) {
          this.httpError = data.message || "Run in shell is not available for this connection.";
          return;
        }
        if (r.status === 400 && data.code === "VALIDATION_FAILED") {
          this.validateResult = data.validation;
          this.httpError = "Validation failed.";
          return;
        }
        if (!r.ok) {
          this.httpError = data.message || r.statusText;
          return;
        }
        this.launchResult = data;
      } finally {
        this.loading = false;
      }
    },

    confirmQuit() {
      this.httpError = "";
      this.notice = "";
      if (!this.launchSupported) {
        this.httpError =
          "Quit server is only available when Launchpad is opened from the same machine as the server (localhost).";
        return;
      }
      const ok = window.confirm(
        "Quit the STAR Server?\n\nThis will stop Launchpad and MCP on this host."
      );
      if (!ok) return;
      void this.doQuitServer();
    },

    async doQuitServer() {
      this.loading = true;
      this.httpError = "";
      this.notice = "";
      try {
        const api = launchpadApiBase();
        const r = await fetch(`${api}/quit`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ confirm: true }),
        });
        let data = null;
        try {
          data = await r.json();
        } catch {
          data = null;
        }
        if (!r.ok) {
          this.httpError = (data && data.message) ? data.message : r.statusText;
          return;
        }
        this.notice = (data && data.message) ? data.message : "Server is shutting down.";
      } finally {
        this.loading = false;
      }
    },

    handleGlobalKey(e) {
      const key = String(e.key || "").toLowerCase();
      if (e.ctrlKey && e.shiftKey && key === "q") {
        e.preventDefault();
        this.confirmQuit();
      }
    },

    slParseInputValues() {
      const raw = (this.slInputJson || "").trim();
      if (!raw) return null;
      const obj = JSON.parse(raw);
      if (!obj || typeof obj !== "object" || Array.isArray(obj)) {
        throw new Error("input_values JSON must be a non-array object");
      }
      return obj;
    },

    async slPostJson(path, payload) {
      const api = launchpadApiBase();
      const r = await fetch(`${api}${path}`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
      const data = await r.json().catch(() => ({}));
      return { r, data };
    },

    async slTranslate() {
      this.loading = true;
      this.httpError = "";
      this.slInvalidateScriptLaneDerivedState();
      try {
        const { r, data } = await this.slPostJson("/script-lane/translate", {
          script_text: this.slScript,
        });
        if (!r.ok) {
          this.httpError = data.message || r.statusText;
          return;
        }
        this.slTranslateResult = data;
        if (data && data.ok && data.ir) {
          this.slIr = data.ir;
          this.slNormalizedScript = data.normalized_script || null;
        }
      } catch (e) {
        this.httpError = e && e.message ? e.message : String(e);
      } finally {
        this.loading = false;
      }
    },

    async slViability() {
      if (!this.slIr) {
        this.httpError = "Translate successfully first so an IR is available.";
        return;
      }
      this.loading = true;
      this.httpError = "";
      this.slViabilityResult = null;
      this.slExecuteResult = null;
      try {
        let input_values = undefined;
        try {
          const parsed = this.slParseInputValues();
          if (parsed && Object.keys(parsed).length) input_values = parsed;
        } catch (e) {
          this.httpError = e && e.message ? e.message : String(e);
          return;
        }
        const body = { ir: this.slIr };
        if (input_values !== undefined) body.input_values = input_values;
        const { r, data } = await this.slPostJson("/script-lane/viability", body);
        if (!r.ok) {
          this.httpError = data.message || r.statusText;
          return;
        }
        this.slViabilityResult = data;
      } catch (e) {
        this.httpError = e && e.message ? e.message : String(e);
      } finally {
        this.loading = false;
      }
    },

    async slExecute() {
      if (!this.slIr) {
        this.httpError = "Translate successfully first so an IR is available.";
        return;
      }
      const wd = (this.slWorkdir || "").trim();
      if (!wd) {
        this.httpError = "Set workdir to an existing or creatable directory on the server host.";
        return;
      }
      let input_values;
      try {
        input_values = this.slParseInputValues();
      } catch (e) {
        this.httpError = e && e.message ? e.message : String(e);
        return;
      }
      if (!input_values || !Object.keys(input_values).length) {
        this.httpError = "input_values JSON must be a non-empty object for execute.";
        return;
      }
      this.loading = true;
      this.httpError = "";
      this.slExecuteResult = null;
      try {
        const { r, data } = await this.slPostJson("/script-lane/execute", {
          ir: this.slIr,
          input_values,
          workdir: wd,
        });
        if (!r.ok) {
          this.httpError = data.message || r.statusText;
          return;
        }
        this.slExecuteResult = data;
      } catch (e) {
        this.httpError = e && e.message ? e.message : String(e);
      } finally {
        this.loading = false;
      }
    },

    slIrJsonText() {
      if (!this.slIr) return "";
      return JSON.stringify(this.slIr, null, 2);
    },

    async slCopyIr() {
      const t = this.slIrJsonText();
      if (!t) return;
      try {
        await navigator.clipboard.writeText(t);
      } catch (e) {
        this.httpError = "Clipboard unavailable: " + (e && e.message ? e.message : String(e));
      }
    },

    slDownloadIr() {
      const t = this.slIrJsonText();
      if (!t) return;
      const blob = new Blob([t], { type: "application/json" });
      const a = document.createElement("a");
      const url = URL.createObjectURL(blob);
      a.href = url;
      a.download = "script-lane-ir.json";
      a.click();
      URL.revokeObjectURL(url);
    },

    slTemporalFallbackBanner() {
      const tr = this.slTranslateResult;
      if (tr && tr.ok === false && tr.fallback_to_temporal) return true;
      const v = this.slViabilityResult;
      if (v && v.ok === false && v.fallback_to_temporal) return true;
      const x = this.slExecuteResult;
      if (x && x.ok === false && x.fallback_to_temporal) return true;
      return false;
    },
  };
}
