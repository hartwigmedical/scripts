#!/usr/bin/env zsh

set -o pipefail

CONFIG_DIR="$(readlink -f "$(dirname "$0")")/.tunnel_configurations"

function print_available_configs() {
    echo "Known configurations:"
    for c in "$(ls $CONFIG_DIR)"; do
        echo "${c}" | sed -e 's/^/  /' -e 's/\.config$//'
    done
}


function print_usage() {
    cat<<EOM

Establish a connection to an application running in a private Kubernetes cluster    
USAGE: $0 (application) [-c config_dir] [-n namespace] [-p local_port]

  application - Application to connect to, see below for existing list
  config_dir  - [optional] directory to use for configuration files instead of the default [.../$(basename $CONFIG_DIR)]
  namespace   - [optional] namespace overriding any set in configuration or the default value "default"
  local_port  = [optional] local port to bind, used to avoid conflicts (overrides default in configuragion file)

NOTE: To start multiple tunnels  use `-p` to avoid local port conflicts and multiple terminals or Ctrl-Z to run many instances

EOM

  print_available_configs
}

gcv="$(gcloud --version | head -n1 | awk '{print $NF}' | awk -F. '{print $1}')"
[[ $gcv -lt "496" ]] && echo "Local version of gcloud is not new enough (need >= 496.0.0)" && exit 1
[[ $# -eq 0 ]] && echo "ERROR: Application must be provided!" && print_usage && exit 1
app="$1" && shift
config_dir="$CONFIG_DIR"
while getopts 'c:n:p:' opt; do
    case "$opt" in
        c)
            config_dir="$OPTARG"
            ;;
        n)
            user_ns="$OPTARG"
            ;;
        p)
            local_port="$OPTARG"
            ;;
        *)
            echo "Unknown argument [$opt]"
            print_usage
            exit 1
            ;;
    esac
done

shift "$(($OPTIND -1))"

config_file="${config_dir}/${app}.config"
if [[ ! -f $config_file ]]; then
    echo "Unknown app [$app]: [$config_file] does not exist"
    print_available_configs 
    exit 1
fi

source $config_file
ns="${user_ns:-$NAMESPACE}"
ns="${ns:-default}"

for var in PROJECT PORT CONTEXT_ROOT CLUSTER; do
    [[ -z ${(P)var} ]] && echo "$var must be specified in [$config_file]" && exit 1
done

local_port="${local_port:-$PORT}"
pod_name_prefix="${POD:=actin-$app}"

gcloud container clusters get-credentials "$CLUSTER" --zone europe-west4 --project "$PROJECT" --dns-endpoint

pod="$(kubectl -n "$ns" get pods | grep "^${pod_name_prefix}-" | awk '{print $1}' | tail -n1)"
[[ -z "$pod" ]] && echo "No pod found for application '$app', check [$config_file]" && exit 1
username="$(gcloud compute os-login describe-profile --format="json" | jq -r '.posixAccounts[].username')"
(sleep 3; open "http://localhost:${local_port}/${CONTEXT_ROOT}?username=${username}" ) &
echo "Establishing tunnel, Ctrl-C will close it"
kubectl -n "$ns" port-forward "pod/${pod}" "$local_port:$PORT"

