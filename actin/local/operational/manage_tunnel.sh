#!/usr/bin/env bash

REMOTE_TUNNEL_PORT=8888
CONFIG_DIR="$(readlink -f "$(dirname "$0")")/.tunnel_configurations"

function print_available_configs() {
  echo "Known configurations:"
  for c in "$(ls $CONFIG_DIR)"; do
    echo "${c}" | sed -e 's/^/  /' -e 's/\.config$//'
  done
}

function cleanup() {
    port=$1
    echo "Cleaning up old processes listening on port $port"
    ps aux | grep "$port:localhost:$port" | grep 'compute start-iap-tunnel' | awk '{print $2}' | while read pid; do
        [[ -n $pid ]] && kill $pid 
        echo "  Process [$pid] killed"
    done
    echo "Cleanup complete"
}

function k8() {
    port=$1 && shift 
    HTTPS_PROXY=localhost:$port kubectl $@
}

if [[ $# -lt 2 ]]; then
    cat<<EOM
    
USAGE: $0 (config name) start|stop [local port]
  config name - Configuration file to read, see below for existing list
  start|stop  - Action to perform
  local port  - Optional local port number to use for the tunnel to avoid conflicts

EOM
print_available_configs
    exit 1
fi

config_file="${CONFIG_DIR}/${1}.config"
[[ ! -f $config_file ]] && echo "Unknown configuration \"$1\": no file [$(basename ${config_file})] in [${CONFIG_DIR}]" && print_available_configs && exit 1

source $config_file

for var in PROJECT CLUSTER POD_NAME REMOTE_APPLICATION_PORT LOCAL_APPLICATION_PORT; do
  [[ -z ${!var} ]] && echo "$var must be specified in [$config_file]" && exit 1
done

if [[ "$PROJECT" == "actin-emc" ]]; then
  [[ -z ${NAMESPACE} ]] && echo "NAMESPACE must be specified in [$config_file] when PROJECT is actin-emc" && exit 1
fi

namespace=""
[[ -n "$NAMESPACE" ]] && namespace="--namespace $NAMESPACE"

action="$2"
[[ $action != "start" && $action != "stop" ]] && echo "Provide a verb (either \"start\" or \"stop\")" && exit 1

GC="gcloud --project $PROJECT compute"

local_tunnel_port="${3:-$REMOTE_TUNNEL_PORT}"

cleanup $local_tunnel_port
if [[ $action == "stop" ]]; then
    echo "Tunnel stopped"
else
    echo "Starting tunnel on port $local_tunnel_port"
    echo "Fetching cluster credentials..."
    gcloud container clusters get-credentials "$CLUSTER" --zone europe-west4 --project "$PROJECT"
    echo "Attempting to determine name of bastion VM"
    read instance zone <<< $($GC instances list | grep bastion | head -n1 | awk '{print $1, $2}')
    [[ -z $instance || -z $zone ]] && echo "Cannot locate bastion VM instance" && exit 1
    echo "Establishing tunnel to $instance"
    $GC ssh $instance --tunnel-through-iap --zone $zone -- -L "$local_tunnel_port:localhost:$REMOTE_TUNNEL_PORT" -N -q -f
    echo "Forwarding $REMOTE_APPLICATION_PORT to localhost:$LOCAL_APPLICATION_PORT for ${POD_NAME}"
    k8 $local_tunnel_port port-forward $(k8 $port get pods $namespace | grep $POD_NAME | awk '{print $1}') \
        "$LOCAL_APPLICATION_PORT:$REMOTE_APPLICATION_PORT" &
    echo "Tunnel started; access $POD_NAME at http://localhost:$LOCAL_APPLICATION_PORT"
fi
