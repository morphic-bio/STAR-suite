const vscode = require('vscode');
const path = require('path');

const EXTENSION_ID = 'morphic-bio.star-suite-walkthrough';
const ROOT_WALKTHROUGH_ID = `${EXTENSION_ID}#setupReference`;

function workspaceRoot() {
  const folder = vscode.workspace.workspaceFolders && vscode.workspace.workspaceFolders[0];
  return folder ? folder.uri.fsPath : null;
}

async function openRelative(relPath) {
  const root = workspaceRoot();
  if (!root) {
    vscode.window.showErrorMessage('STAR-suite walkthrough requires an open workspace.');
    return;
  }
  const target = vscode.Uri.file(path.join(root, relPath));
  const doc = await vscode.workspace.openTextDocument(target);
  await vscode.window.showTextDocument(doc, { preview: false });
}

function launchTerminal(name, command) {
  const root = workspaceRoot();
  if (!root) {
    vscode.window.showErrorMessage('STAR-suite walkthrough requires an open workspace.');
    return;
  }
  const terminal = vscode.window.createTerminal({ name, cwd: root });
  terminal.show(true);
  terminal.sendText(command, true);
}

function quote(arg) {
  return `'${String(arg).replace(/'/g, `'"'"'`)}'`;
}

function bashFromRepo(command) {
  const root = workspaceRoot();
  if (!root) {
    return null;
  }
  return `bash -lc ${quote(`cd ${quote(root)} && ${command}`)}`;
}

function register(context, command, fn) {
  context.subscriptions.push(vscode.commands.registerCommand(command, fn));
}

function activate(context) {
  register(context, 'starSuiteWalkthrough.openMainWalkthrough', async () => {
    try {
      await vscode.commands.executeCommand('workbench.action.openWalkthrough', ROOT_WALKTHROUGH_ID, true);
    } catch (error) {
      const detail = error instanceof Error ? error.message : String(error);
      vscode.window.showInformationMessage(`Open Getting Started and select STAR-suite Demo Walkthrough. ${detail}`);
    }
  });

  register(context, 'starSuiteWalkthrough.openOverviewGuide', () => openRelative('docs/codespaces/00_overview.md'));
  register(context, 'starSuiteWalkthrough.openSetupGuide', () => openRelative('docs/codespaces/01_setup_reference.md'));
  register(context, 'starSuiteWalkthrough.openBulkGuide', () => openRelative('docs/codespaces/02_bulk.md'));
  register(context, 'starSuiteWalkthrough.openSlamGuide', () => openRelative('docs/codespaces/03_slam.md'));
  register(context, 'starSuiteWalkthrough.openFixtureGuide', () => openRelative('docs/codespaces/04_single_cell_fixture.md'));
  register(context, 'starSuiteWalkthrough.openPerturbGuide', () => openRelative('docs/codespaces/05_perturb.md'));
  register(context, 'starSuiteWalkthrough.openFlexGuide', () => openRelative('docs/codespaces/06_flex.md'));

  register(context, 'starSuiteWalkthrough.runFetchReference', () => {
    const cmd = bashFromRepo('bash scripts/codespaces/fetch_public_chr22y_reference.sh');
    if (cmd) {
      launchTerminal('STAR-suite Walkthrough: Reference Fetch', cmd);
    }
  });

  register(context, 'starSuiteWalkthrough.runBuildIndex', () => {
    const cmd = bashFromRepo('bash scripts/codespaces/build_public_chr22y_index.sh --threads 4');
    if (cmd) {
      launchTerminal('STAR-suite Walkthrough: Build Index', cmd);
    }
  });

  register(context, 'starSuiteWalkthrough.runFixtureDryRun', () => {
    const cmd = bashFromRepo('bash scripts/codespaces/derive_public_chr22y_gex_fixture.sh --threads 4 --read-limit 200000 --dry-run');
    if (cmd) {
      launchTerminal('STAR-suite Walkthrough: Fixture Dry Run', cmd);
    }
  });

  register(context, 'starSuiteWalkthrough.runBulkDryRun', () => {
    const cmd = bashFromRepo('bash scripts/codespaces/run_bulk_public_demo.sh --dry-run --genome-dir /path/to/star_bulk_index');
    if (cmd) {
      launchTerminal('STAR-suite Walkthrough: Bulk Dry Run', cmd);
    }
  });

  register(context, 'starSuiteWalkthrough.runSlamDryRun', () => {
    const cmd = bashFromRepo('bash scripts/codespaces/run_slam_public_demo.sh --dry-run --genome-dir /path/to/star_bulk_index');
    if (cmd) {
      launchTerminal('STAR-suite Walkthrough: SLAM Dry Run', cmd);
    }
  });

  register(context, 'starSuiteWalkthrough.runPerturbDryRun', () => {
    const cmd = bashFromRepo('bash scripts/codespaces/run_perturb_public_demo.sh --dry-run');
    if (cmd) {
      launchTerminal('STAR-suite Walkthrough: Perturb Dry Run', cmd);
    }
  });

  register(context, 'starSuiteWalkthrough.runFlexDryRun', () => {
    const cmd = bashFromRepo('bash scripts/codespaces/run_flex_public_demo.sh --dry-run');
    if (cmd) {
      launchTerminal('STAR-suite Walkthrough: Flex Dry Run', cmd);
    }
  });
}

function deactivate() {}

module.exports = {
  activate,
  deactivate,
};
