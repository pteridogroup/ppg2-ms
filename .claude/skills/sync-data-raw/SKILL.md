---
name: sync-data-raw
description: Bidirectionally sync the ppg2-ms project's data_raw/ folder (raw data, not under version control, containing personally identifiable information) between Joel's Mac and the nittalab server, keeping the newest version of each file by content. Use when Joel asks to sync data_raw, mentions data_raw has drifted between machines, or wants to reconcile raw data files between his Mac and nittalab.
user-invocable: true
---

# sync-data-raw — reconcile data_raw/ between Mac and nittalab

`data_raw/` holds raw data for the "pre" workflow (`_targets_pre.R`) — things like
contributor emails and voting spreadsheets that include personally identifiable
information, so they're deliberately **not** committed to git (see `.gitignore` /
`.dockerignore`) or synced via GitHub. It only exists in two places, and they can drift
out of sync since each is edited independently:

- **Mac**: `/Users/joelnitta/repos/ppg2-ms/data_raw` — a symlink into Dropbox, resolving
  to `/Users/joelnitta/Library/CloudStorage/Dropbox/project_data/ppg2_ms/data_raw`
- **nittalab**: `/home/jnitta/ppg2-ms/data_raw` — a plain directory (**not** the
  devcontainer git clone at `~/repos/ppg2-ms`, which is a different checkout entirely)

## What to exclude

- `.DS_Store` — macOS Finder junk, not real data
- `ppg-0.0.0.9003/` — a cached snapshot of the `ppg` R package source tree; low value to
  fuss over syncing since it can just be re-fetched from GitHub if ever needed

Both are excluded by default below. If the user adds genuinely new raw-data files later,
extend the `--exclude` list rather than removing these two.

## Workflow

1. **Resolve the Mac-side path** (the symlink target, since rsync needs a real path):
   ```
   MAC_DIR="$(realpath /Users/joelnitta/repos/ppg2-ms/data_raw)/"
   ```

2. **Dry-run both directions first** and show the user what would change before touching
   anything — `--itemize-changes` output distinguishes a real content change (flagged
   `c`/`s` for checksum/size differing) from a harmless metadata-only touch (mtime/perm
   sync on already-identical content):
   ```
   rsync -au --checksum --exclude='.DS_Store' --exclude='ppg-0.0.0.9003/' \
     --itemize-changes --dry-run -e ssh "$MAC_DIR" nittalab:/home/jnitta/ppg2-ms/data_raw/

   rsync -au --checksum --exclude='.DS_Store' --exclude='ppg-0.0.0.9003/' \
     --itemize-changes --dry-run -e ssh nittalab:/home/jnitta/ppg2-ms/data_raw/ "$MAC_DIR"
   ```
   `-u` (update) skips a file if the destination's mtime is newer; `--checksum` forces a
   real content comparison rather than trusting size+mtime, so files that only *look*
   different (e.g. re-copied at different times but byte-identical) aren't needlessly
   retransferred, and files that genuinely differ get caught even if their size happens
   to match.

3. **Confirm with the user** if anything beyond a metadata-only touch shows up —
   especially if the *same* file appears changed in *both* directions (a real conflict:
   both copies were edited since they last matched, and one edit will be silently
   discarded). That shouldn't happen under normal one-writer-at-a-time use, but check for
   it rather than assuming.

4. **Run for real** (drop `--dry-run`) once confirmed:
   ```
   rsync -au --checksum --exclude='.DS_Store' --exclude='ppg-0.0.0.9003/' \
     --itemize-changes -e ssh "$MAC_DIR" nittalab:/home/jnitta/ppg2-ms/data_raw/

   rsync -au --checksum --exclude='.DS_Store' --exclude='ppg-0.0.0.9003/' \
     --itemize-changes -e ssh nittalab:/home/jnitta/ppg2-ms/data_raw/ "$MAC_DIR"
   ```

5. **Verify** by checksumming the top-level files on both sides and confirming they
   match:
   ```
   shasum -a 256 "$MAC_DIR"*.csv
   ssh nittalab 'sha256sum /home/jnitta/ppg2-ms/data_raw/*.csv'
   ```

## Notes

- This is a plain two-way "newest content wins" merge — it never deletes files. If a file
  was deliberately removed on one side, this will resurrect it from the other side rather
  than propagating the deletion. Flag that to the user if you notice a file present on
  only one side (could be a genuine deletion, or just a new file not yet synced) rather
  than assuming either way.
- Requires SSH access to `nittalab` (see `~/.ssh/config`); if the connection times out,
  it's usually a VPN issue on Joel's end (private `10.19.x.x` address) — ask him to check
  before assuming something is broken.
