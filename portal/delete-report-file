#!/usr/bin/env bash

set -euo pipefail

# Delete a manually added report file from a report in REVIEW status.
#
# Optional environment variables:
# - SERVER_URL: Portal base URL, for example http://portal.gateway.pilot-1/ or http://localhost:8082
# - AUTH_TOKEN: Bearer token from the auth service offline-session flow
# - AUTH_URL: Override the offline-session token URL
# - SUBMISSION_ID: Portal submission id (preferred)
# - REPORT_ID: Portal report id (used when SUBMISSION_ID is omitted)
# - COMPONENT_ID: Component id to delete. If omitted, the script lists components first.

prompt() {
  local var_name="$1"
  local label="$2"
  local default_value="${3:-}"
  local value

  if [ -n "${!var_name:-}" ]; then
    return
  fi

  if [ -n "$default_value" ]; then
    read -r -p "$label [$default_value]: " value
    printf -v "$var_name" '%s' "${value:-$default_value}"
  else
    read -r -p "$label: " value
    printf -v "$var_name" '%s' "$value"
  fi
}

extract_token() {
  python3 -c '
from urllib.parse import parse_qs
import sys
s = sys.argv[1]
print(parse_qs(s.partition("#")[2]).get("refreshToken", [s])[0])
' "$1"
}

print_components() {
  python3 -c '
import json
import sys

components = json.load(sys.stdin)

if not components:
    print("No components found.")
    sys.exit(0)

print("{:>8}  {:<8}  {:>12}  {}".format("ID", "TYPE", "SIZE", "PATH"))

for component in components:
    folder = component.get("folder") or ""
    name = component.get("name") or ""
    path = f"{folder}/{name}" if folder else name

    component_id = component.get("id", "")
    file_type = component.get("fileType", "")
    size = component.get("sizeInBytes", "")

    print("{:>8}  {:<8}  {:>12}  {}".format(
        component_id,
        file_type,
        size,
        path,
    ))
'
}

prompt SERVER_URL "Portal server URL" "http://portal.gateway.pilot-1/"
SERVER_URL="${SERVER_URL%/}"

REPORT_OR_SUBMISSION_ID="${SUBMISSION_ID:-${REPORT_ID:-}}"
prompt REPORT_OR_SUBMISSION_ID "Submission id or report id"

if [ -z "${AUTH_TOKEN:-}" ]; then
  AUTH_URL="${AUTH_URL:-$SERVER_URL/auth/offline-session/token?rd=/api/roles}"
  echo
  echo "Open this URL, authenticate, and copy the token from the redirect URL fragment:"
  echo "$AUTH_URL"
  echo
  read -r -p "Paste token or full redirect URL: " pasted_token
  AUTH_TOKEN="$(extract_token "$pasted_token")"
fi

if [ -z "$AUTH_TOKEN" ]; then
  echo "Error: AUTH_TOKEN is empty" >&2
  exit 1
fi

if [[ "$REPORT_OR_SUBMISSION_ID" == *-* ]]; then
  echo "Resolving report id for submission $REPORT_OR_SUBMISSION_ID..."
  result_response="$(
    curl -fsS \
      -H "Authorization: Bearer $AUTH_TOKEN" \
      "$SERVER_URL/api/submission/$REPORT_OR_SUBMISSION_ID/result"
  )"
  REPORT_ID="$(
    python3 -c 'import json, sys; print(json.load(sys.stdin)["reportId"])' <<< "$result_response"
  )"
else
  REPORT_ID="$REPORT_OR_SUBMISSION_ID"
fi

if ! [[ "$REPORT_ID" =~ ^[0-9]+$ ]]; then
  echo "Error: resolved report id must be numeric: $REPORT_ID" >&2
  exit 1
fi

if [ -z "${COMPONENT_ID:-}" ]; then
  echo "Fetching components for report $REPORT_ID..."
  components="$(
    curl -fsS \
      -H "Authorization: Bearer $AUTH_TOKEN" \
      "$SERVER_URL/api/reports/$REPORT_ID/components"
  )"
  printf '%s\n' "$components" | print_components
  echo
  prompt COMPONENT_ID "Component id to delete"
fi

if ! [[ "$COMPONENT_ID" =~ ^[0-9]+$ ]]; then
  echo "Error: component id must be numeric: $COMPONENT_ID" >&2
  exit 1
fi

read -r -p "Delete component $COMPONENT_ID from report $REPORT_ID? [y/N]: " confirm
case "$confirm" in
  y|Y|yes|YES) ;;
  *)
    echo "Cancelled."
    exit 0
    ;;
esac

echo "Deleting component $COMPONENT_ID from report $REPORT_ID..."
curl -fsS -X DELETE \
  -H "Authorization: Bearer $AUTH_TOKEN" \
  "$SERVER_URL/api/reports/$REPORT_ID/components/$COMPONENT_ID"

echo
echo "Done: deleted component $COMPONENT_ID from report $REPORT_ID."
