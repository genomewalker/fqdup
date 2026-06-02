#!/usr/bin/env bash
# Usage: scripts/release.sh <patch|minor|major|x.y.z> [-y]
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION_FILE="$REPO_ROOT/include/fqdup/version.hpp"
CHANGELOG="$REPO_ROOT/CHANGELOG.md"

BUMP="${1:-}"
AUTO_YES=0
[[ "${2:-}" == "-y" ]] && AUTO_YES=1

if [[ -z "$BUMP" ]]; then
    echo "Usage: $0 <patch|minor|major|x.y.z> [-y]" >&2
    exit 1
fi

CURRENT=$(grep -oP '(?<=FQDUP_VERSION ")[\d.]+(?=")' "$VERSION_FILE")
IFS='.' read -r MAJ MIN PAT <<< "$CURRENT"

case "$BUMP" in
    major) NEW="$((MAJ+1)).0.0" ;;
    minor) NEW="${MAJ}.$((MIN+1)).0" ;;
    patch) NEW="${MAJ}.${MIN}.$((PAT+1))" ;;
    [0-9]*.[0-9]*.[0-9]*) NEW="$BUMP" ;;
    *) echo "Invalid bump: $BUMP" >&2; exit 1 ;;
esac

echo "Release: $CURRENT → $NEW"

if [[ $AUTO_YES -eq 0 ]]; then
    read -r -p "Proceed? [y/N] " REPLY
    [[ "$REPLY" =~ ^[Yy]$ ]] || exit 0
fi

cd "$REPO_ROOT"

if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "Error: uncommitted changes present" >&2; exit 1
fi

if git rev-parse "v$NEW" &>/dev/null; then
    echo "Error: tag v$NEW already exists" >&2; exit 1
fi

TODAY=$(date +%Y-%m-%d)

# Bump version header
sed -i "s/FQDUP_VERSION \"$CURRENT\"/FQDUP_VERSION \"$NEW\"/" "$VERSION_FILE"

# Update or prepend CHANGELOG entry
if grep -q "^## \[$NEW\]" "$CHANGELOG"; then
    # Entry exists — stamp today's date
    sed -i "s|^## \[$NEW\].*|## [$NEW] - $TODAY|" "$CHANGELOG"
else
    # Insert blank entry after the first line
    TMP=$(mktemp)
    head -1 "$CHANGELOG" > "$TMP"
    printf '\n## [%s] - %s\n\n_Release notes: see commits since v%s._\n' \
        "$NEW" "$TODAY" "$CURRENT" >> "$TMP"
    tail -n +2 "$CHANGELOG" >> "$TMP"
    mv "$TMP" "$CHANGELOG"
fi

git add "$VERSION_FILE" "$CHANGELOG"
git commit -m "release: v$NEW"
git tag "v$NEW"
git push origin main
git push origin "v$NEW"

echo "Released v$NEW — CI will create the GitHub release automatically."
