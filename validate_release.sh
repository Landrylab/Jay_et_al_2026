#!/usr/bin/env bash
# =============================================================================
# Description
#   Validate shell, Python, and R syntax across the shareable bundle.
# =============================================================================

set -euo pipefail

script_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
status=0

while IFS= read -r path; do
  bash -n "${path}" || status=1
done < <(find "${script_root}" -type f -name '*.sh' -print)

while IFS= read -r path; do
  python3 -c 'import ast, pathlib, sys; ast.parse(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))' "${path}" || status=1
done < <(find "${script_root}" -type f -name '*.py' -print)

while IFS= read -r path; do
  Rscript -e 'parse(file=commandArgs(TRUE)[1])' "${path}" >/dev/null || status=1
done < <(find "${script_root}" -type f -name '*.R' -print)

scan_options=(
  --recursive
  --line-number
  --ignore-case
  --extended-regexp
  --include='*.sh'
  --include='*.py'
  --include='*.R'
  --include='*.md'
  --exclude='validate_release.sh'
)

if grep "${scan_options[@]}" \
  '/Users/|/home/[^/]+/|[A-Za-z]:\\Users\\|Downloads/|Desktop/' \
  "${script_root}"; then
  echo "Personal or workstation-specific path found in the public bundle." >&2
  status=1
fi

if grep "${scan_options[@]}" \
  '(^|[^[:alnum:]_])(you|your|yours|yourself|yourselves)($|[^[:alnum:]_])' \
  "${script_root}"; then
  echo "Second-person wording found in the public bundle." >&2
  status=1
fi

if grep "${scan_options[@]}" \
  'as an AI|AI-generated|delve|leverage|utilize' \
  "${script_root}"; then
  echo "Prohibited generic or AI-style wording found in the public bundle." >&2
  status=1
fi

exit "${status}"
