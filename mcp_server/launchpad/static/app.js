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

function launchpadApp() {
  return {
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
    httpError: "",
    notice: "",
    showAbout: false,

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

    paramHasRichTip(pname) {
      const p = this.paramByName(pname);
      if (!p) return false;
      if (String(p.description || "").trim()) return true;
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
      } else {
        this.launchSupported = false;
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

    async doValidate() {
      this.loading = true;
      this.httpError = "";
      try {
        const api = launchpadApiBase();
        const r = await fetch(
          `${api}/workflows/${encodeURIComponent(this.workflowId)}/validate`,
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
        this.validateResult = data;
      } finally {
        this.loading = false;
      }
    },

    async doRender() {
      await this.doValidate();
      if (!this.validateResult?.valid) return;
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
      await this.doValidate();
      if (!this.validateResult?.valid) return;
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
            body: JSON.stringify({ params: this.collectParams() }),
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
  };
}
