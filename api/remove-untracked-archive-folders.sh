#!/usr/bin/env bash
set -euo pipefail

REPORT_URL="${REPORT_URL:-http://api.pilot-1/hmf/v2/archive/deletion-report}"
BATCH_SIZE="${BATCH_SIZE:-50}"
PARALLELISM="${PARALLELISM:-4}"

usage() {
  cat <<EOF
Usage: $(basename "$0")

Fetches the hartwig-api archive deletion report, prompts for confirmation,
then deletes only folders listed in .untrackedFolders.

Environment variables:
  REPORT_URL   Deletion report URL. Default: $REPORT_URL
  BATCH_SIZE   Number of folder URIs per gcloud invocation. Default: $BATCH_SIZE
  PARALLELISM  Number of parallel gcloud invocations. Default: $PARALLELISM

Example:
  REPORT_URL=http://api.pilot-1/hmf/v2/archive/deletion-report \\
  BATCH_SIZE=50 \\
  PARALLELISM=4 \\
  $0
EOF
}

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if ! [[ "$BATCH_SIZE" =~ ^[1-9][0-9]*$ ]]; then
  echo "BATCH_SIZE must be a positive integer, got: $BATCH_SIZE" >&2
  exit 1
fi

if ! [[ "$PARALLELISM" =~ ^[1-9][0-9]*$ ]]; then
  echo "PARALLELISM must be a positive integer, got: $PARALLELISM" >&2
  exit 1
fi

require_command curl
require_command jq
require_command gcloud

work_dir="$(mktemp -d "${TMPDIR:-/tmp}/archive-report-cleanup.XXXXXX")"
report_file="$work_dir/deletion-report.json"
targets_file="$work_dir/untracked-folder-uris.txt"

cleanup() {
  rm -rf "$work_dir"
}
trap cleanup EXIT

echo "Fetching archive deletion report from: $REPORT_URL"
curl -fsS "$REPORT_URL" -o "$report_file"

jq -e '
  has("status")
  and has("scannedFolders")
  and has("trackedFolders")
  and (.untrackedFolders | type == "array")
  and (.staleDeletions | type == "array")
  and (.expiredPendingDeletionTasks | type == "array")
' "$report_file" >/dev/null

jq -r '
  .untrackedFolders[]
  | "\(.archiveUri)/\(.submissionId)/\(.analysisId)/\(.category)"
' "$report_file" > "$targets_file"

target_count="$(wc -l < "$targets_file" | tr -d ' ')"

echo
jq -r '
  "Report status: \(.status)",
  "Scanned folders: \(.scannedFolders)",
  "Tracked folders: \(.trackedFolders)",
  "Untracked folders: \(.untrackedFolders | length)",
  "Stale deletions: \(.staleDeletions | length)",
  "Expired pending deletion tasks: \(.expiredPendingDeletionTasks | length)"
' "$report_file"

if [[ "$target_count" -eq 0 ]]; then
  echo
  echo "No untracked folders reported; nothing to delete."
  exit 0
fi

echo
echo "Untracked folders by archive URI:"
jq -r '.untrackedFolders[].archiveUri' "$report_file" | sort | uniq -c

echo
echo "Untracked folders by category:"
jq -r '.untrackedFolders[].category' "$report_file" | sort | uniq -c

echo
echo "Full target list:"
sed 's/^/  /' "$targets_file"

echo
echo "About to delete $target_count untracked archive folder(s)."
echo "Batch size: $BATCH_SIZE; parallel gcloud invocations: $PARALLELISM"
read -r -p "Type DELETE to continue: " confirmation

if [[ "$confirmation" != "DELETE" ]]; then
  echo "Aborted; no folders were deleted."
  exit 1
fi

echo "Deleting reported untracked archive folders..."
xargs -n "$BATCH_SIZE" -P "$PARALLELISM" \
  gcloud storage rm --recursive --continue-on-error \
  < "$targets_file"

echo "Finished deleting reported untracked archive folders."
