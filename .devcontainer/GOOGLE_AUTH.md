# Google Sheets authentication (googlesheets4)

`_targets_pre.R` reads several Google Sheets (`taxon_comments`, votes,
emails) via `googlesheets4::read_sheet()`. This needs a cached OAuth
token on disk — there's no service account, so `gs4_auth()` must be
run interactively at least once per machine.

## Symptom

```
Error in tar_make():
  Can't get Google credentials.
```

This happens when there's no cached token in `~/.cache/gargle` — that
is gargle's actual default (`rappdirs::user_cache_dir("gargle")`),
*not* `~/.cache/R/gargle`. Inside the devcontainer, that directory
used to live only in the container's ephemeral home dir, so it was
wiped on every rebuild — hence needing to redo this after each
rebuild.

## Fix (one-time per machine, persists across rebuilds)

Both `.devcontainer/devcontainer.json` and `.devcontainer/mac/devcontainer.json`
bind-mount `~/.cache/gargle` from the host, the same way they already
mount `~/.Renviron`. As long as the directory exists on the host, the
cached token survives container rebuilds.

1. **On the host** (lab server or Mac, outside the container), create
   the directory once so Docker doesn't auto-create it as root-owned:
   ```sh
   mkdir -p ~/.cache/gargle
   ```
2. **Rebuild the devcontainer** (VS Code: "Dev Containers: Rebuild
   Container") so it picks up the mount.
3. **Inside the container**, authenticate interactively — not via
   `tar_make()`, which runs in a non-interactive subprocess and can't
   open a browser:
   ```r
   googlesheets4::gs4_auth()
   ```
   This opens a browser via the `$BROWSER` redirect set up in
   `.Rprofile`, and caches the token into `~/.cache/gargle`.
4. Run `tar_make()` — it will reuse the cached token from now on,
   including after future rebuilds.

## Running outside the devcontainer (e.g. laptop, no Docker)

The home directory persists normally there, so `googlesheets4::gs4_auth()`
only needs to be run once, ever, per machine.
