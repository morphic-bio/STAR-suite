# Attribution-history correction — 2026-09-13

The owner requested removal of Claude/Anthropic co-author trailers that contradicted the repository's authorship policy. The correction removed 14 commit-message trailers and changed 516 of 835 commit IDs reachable in the GitHub snapshot. It preserved every commit's file tree, human author and committer, timestamps, message content outside the removed trailer/terminal blank lines, and parent ordering. Parent IDs changed where necessary. Unaffected commits retain their exact object IDs and signatures; signatures invalidated by a rewritten commit were removed (6 commits).

The operation used a selective Git-object rewrite, verified against every original commit. This avoids a full export/import stripping signatures from unaffected upstream history. No commits were pruned, and no source code was changed by the history rewrite. STAR Suite separately applied the already prepared repository-hygiene commit; that cleanup and these documentation/hooks are newer than the historical release tags.

## Historical provenance and releases

- [COMMIT_ID_MAP_20260913.tsv](COMMIT_ID_MAP_20260913.tsv) translates old commit IDs to their rewritten equivalents. Identity rows cover unaffected commits. Search using a full SHA or a uniquely identifying prefix; do not alter old benchmark records or embedded source-revision strings.
- [TAG_ID_MAP_20260913.tsv](TAG_ID_MAP_20260913.tsv) records tag-object changes and unchanged checkout trees. Tag names were retained and updated in place.
- Existing compiled release assets and their SHA-256 checksums are unchanged. A binary may therefore report an old source revision; use the commit map to resolve its equivalent current history. Newly generated GitHub source archives reflect the rewritten Git metadata.
- GitHub mirrors, verified Git bundles, all release assets, and metadata were backed up before rewriting. The owner also retains complete local Git-directory and worktree-edit snapshots, plus local-only commit maps.
- GitHub Actions were temporarily disabled during ref publication to avoid tag-triggered release rebuilds, then restored to their previous settings.

## Existing clones

Fetch the rewritten branches and tags into a fresh clone, or translate existing branch tips through the map after backing up local work. Replaying local commits must preserve their edits and merges. Do not merge the old and rewritten histories together. The owner's known local worktrees were migrated without changing their checked-out source trees; the primary checkout can then fast-forward to the new housekeeping commits while retaining local edits.

GitHub-controlled closed-pull-request refs and cached old commit pages are outside an ordinary force-push. They can retain historical copies even after every writable branch and tag has been corrected; see [GitHub's history-rewrite documentation](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/removing-sensitive-data-from-a-repository). They are included in the backup and map, but are not advertised as writable refs or as erased history.

## Preventing recurrence

Run `python3 .githooks/install.py` after cloning. This installs hooks in the shared hooks directory, covering linked worktrees as well. The installer preserves unrelated existing hooks instead of overwriting them. The commit-message hook rejects Anthropic co-author trailers; the pre-push hook also checks ancestors, preventing a stale local branch from restoring the removed trailers. Product documentation mentioning Claude is unaffected.
