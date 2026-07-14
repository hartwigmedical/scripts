#!/usr/bin/env bash

set -euo pipefail

# Upload a manually added file to a report in REVIEW status.
#
# Optional environment variables:
# - SERVER_URL: Portal base URL, for example http://portal.gateway.pilot-1 or http://localhost:8082
# - AUTH_TOKEN: Bearer token from the auth service offline-session flow
# - AUTH_URL: Override the offline-session token URL
# - LOCAL_FILE: Local file to upload
# - REPORT_FILE_PATH: Report-relative destination path, for example supplementary/example.pdf
# - REPORT_ID: Portal report id

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

json_string() {
  python3 -c 'import json, sys; print(json.dumps(sys.argv[1]))' "$1"
}

url_encode() {
  python3 -c 'import urllib.parse, sys; print(urllib.parse.quote(sys.argv[1], safe=""))' "$1"
}

extract_token() {
  python3 -c '
from urllib.parse import parse_qs
import sys
s = sys.argv[1]
print(parse_qs(s.partition("#")[2]).get("refreshToken", [s])[0])
' "$1"
}

prompt SERVER_URL "Portal server URL" "http://portal.gateway.pilot-1/"
SERVER_URL="${SERVER_URL%/}"

prompt REPORT_ID "Report id"
prompt LOCAL_FILE "Local file path" "$HOME/Desktop/example.pdf"
prompt REPORT_FILE_PATH "Report-relative file path" "$(basename "$LOCAL_FILE")"

if [ ! -f "$LOCAL_FILE" ]; then
  echo "Error: local file does not exist: $LOCAL_FILE" >&2
  exit 1
fi

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

request_body="{\"filePath\":$(json_string "$REPORT_FILE_PATH")}"

echo "Requesting signed upload URL for report $REPORT_ID path '$REPORT_FILE_PATH'..."
upload_response="$(
  curl -fsS -X POST \
    -H "Authorization: Bearer $AUTH_TOKEN" \
    -H "Content-Type: application/json" \
    -d "$request_body" \
    "$SERVER_URL/api/reports/$REPORT_ID/components/upload-url"
)"

signed_url="$(
  python3 -c 'import json, sys; print(json.load(sys.stdin)["signedUrl"])' <<< "$upload_response"
)"

if [ -z "$signed_url" ]; then
  echo "Error: server did not return a signedUrl" >&2
  echo "$upload_response" >&2
  exit 1
fi

echo "Uploading '$LOCAL_FILE'..."
curl -fsS -X PUT \
  -H "Content-Type: application/pdf" \
  --upload-file "$LOCAL_FILE" \
  "$signed_url"

encoded_path="$(url_encode "$REPORT_FILE_PATH")"

echo "Finalizing upload..."
curl -fsS -X POST \
  -H "Authorization: Bearer $AUTH_TOKEN" \
  "$SERVER_URL/api/reports/$REPORT_ID/components/finalize?filePath=$encoded_path"

echo
echo "Done: uploaded '$LOCAL_FILE' as '$REPORT_FILE_PATH' on report $REPORT_ID."
