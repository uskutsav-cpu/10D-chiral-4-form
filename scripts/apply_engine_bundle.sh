#!/usr/bin/env bash
# Apply this release to a clean existing main branch; never force-push or reset.
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  printf 'Usage: bash apply_engine_bundle.sh BUNDLE [EXISTING_REPOSITORY]\n' >&2
  exit 2
fi
BUNDLE_INPUT="$1"
REPO_INPUT="${2:-$HOME/Downloads/10D-chiral-4-form}"
[[ -f "$BUNDLE_INPUT" ]] || { printf 'Bundle not found: %s\n' "$BUNDLE_INPUT" >&2; exit 2; }
BUNDLE="$(cd "$(dirname "$BUNDLE_INPUT")" && pwd)/$(basename "$BUNDLE_INPUT")"
[[ -d "$REPO_INPUT" ]] || { printf 'Repository folder not found: %s\n' "$REPO_INPUT" >&2; exit 2; }
cd "$REPO_INPUT"
ROOT="$(git rev-parse --show-toplevel)"
[[ "$(pwd -P)" == "$(cd "$ROOT" && pwd -P)" ]] || { echo 'Pass the repository root, not a subfolder.' >&2; exit 2; }
ORIGIN="$(git remote get-url origin)"
case "$ORIGIN" in
  https://github.com/uskutsav-cpu/10D-chiral-4-form|https://github.com/uskutsav-cpu/10D-chiral-4-form.git|git@github.com:uskutsav-cpu/10D-chiral-4-form.git) ;;
  *) printf 'Unexpected origin; refusing to modify this repository: %s\n' "$ORIGIN" >&2; exit 2 ;;
esac
[[ "$(git symbolic-ref --quiet --short HEAD)" == main ]] || { echo 'Switch to main first; no files have been modified.' >&2; exit 2; }
[[ -z "$(git status --porcelain)" ]] || { echo 'Working tree has changes or untracked files. Commit/stash them first; nothing will be overwritten.' >&2; exit 2; }
git bundle verify "$BUNDLE"
# Package refs do not change the working tree; incompatible updates are rejected.
git fetch --no-tags "$BUNDLE" \
  refs/heads/main:refs/remotes/engine-package/main \
  refs/heads/research/paper2-computational-engine:refs/remotes/engine-package/research
RELEASE="$(git rev-parse refs/remotes/engine-package/main)"
RESEARCH="$(git rev-parse refs/remotes/engine-package/research)"
BASE=0dedbebd6090017b67ab9999e3d25a1896264e19
git merge-base --is-ancestor "$BASE" "$RELEASE" || { echo 'Package does not extend the expected foundation.' >&2; exit 2; }
git merge-base --is-ancestor HEAD "$RELEASE" || { echo 'Your main contains newer/divergent commits. Stopping without changing files; do not force-push.' >&2; exit 2; }
git merge --ff-only "$RELEASE"
if ! git show-ref --verify --quiet refs/heads/research/paper2-computational-engine; then
  git branch --no-track research/paper2-computational-engine "$RESEARCH"
fi
printf '\nApplied release locally. No remote push was performed.\n'
git log -4 --oneline
printf '\nNext: create/activate .venv, install .[dev], run make verify, then git push origin main.\n'
