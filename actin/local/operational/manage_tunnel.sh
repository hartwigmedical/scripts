#!/usr/bin/env zsh

function print_usage() {
  echo "$0: establish a connection to an application running in Kubernetes"
  echo "Use one of these forms:"
  echo "  $0 prod [app]"
  echo "  $0 research [app] [namespace]"
}

[[ $# -ne 2 && $# -ne 3 ]] && print_usage && exit 1

env="$1"
app="$2"
ns="${3:-default}"

[[ "$env" != "research" && "$env" != "prod" ]] && print_usage && exit 1
[[ "$env" == "prod" && "$ns" != "default" ]] && print_usage && exit 1
[[ "$env" == "research" && "$ns" == "default" ]] && echo "Warning: namespace '$ns' given for research, could be a problem"

declare -A clusters=( ["research"]="research-cluster-patient-processing" ["prod"]="shared-cluster-services" )
declare -A projects=( ["research"]="actin-research" ["prod"]="actin-shared" )
declare -A apps=( ["trial-ui"]="actin-trial-ui" ["tracker-ui"]="actin-patient-tracker" ["feed-review"]="actin-feed-review" )
declare -A ports=( ["trial-ui"]="8082" ["tracker-ui"]="8085" ["feed-review"]="8080" )

set -e

gcloud container clusters get-credentials "${clusters[$env]}" --zone europe-west4 --project "${projects[$env]}" --dns-endpoint

pod="$(kubectl -n "$ns" get pods | grep "^${apps[$app]}-" | awk '{print $1}' | tail -n1)"
[[ -z "$pod" ]] && echo "No pod found for application '$app'" && exit 1
remote="${ports[$app]}"
increment=0
[[ $env == "research" ]] && increment=1000
local="$(( $remote + $increment ))"
username="$(gcloud compute os-login describe-profile --format="json" | jq -r '.posixAccounts[].username')"
(sleep 3; open "http://localhost:$local/${app}?username=${username}" ) &
kubectl -n "$ns" port-forward "pod/${pod}" "$local:$remote"
