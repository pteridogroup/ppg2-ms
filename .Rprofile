source("renv/activate.R")

# In VS Code devcontainers/remote sessions, $BROWSER points to a helper
# script that opens URLs in the local browser (needed for OAuth flows
# like gs4_auth(), since xdg-open isn't installed in the container).
if (nzchar(Sys.getenv("BROWSER"))) {
  options(browser = Sys.getenv("BROWSER"))
}
