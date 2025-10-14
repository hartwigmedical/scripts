#!/usr/bin/env bash

REMOTE_TUNNEL_PORT=8888
CONFIG_DIR="$(readlink -f "$(dirname "$0")")/.tunnel_configurations"

set -x

function print_available_configs() {
  echo "Known configurations:"
  for c in "$(ls $CONFIG_DIR)"; do
    echo "${c}" | sed -e 's/^/  /' -e 's/\.config$//'
  done
}

function cleanup() {
    port=$1
    echo "Cleaning up old processes listening on port $port"
    ps aux | grep "$port:localhost:$REMOTE_TUNNEL_PORT" | grep 'compute start-iap-tunnel' | awk '{print $2}' | while read pid; do
        [[ -n $pid ]] && kill $pid 
        echo "  Process [$pid] killed"
    done
    echo "Cleanup complete"
}

function k8() {
    port=$1 && shift 
    HTTPS_PROXY=localhost:$port kubectl $@
}

function usage() {
    cat<<EOM
    
USAGE: $0 stop [local port] 
       $0 start (config name) [local port]

  start|stop  - Action to perform
  config name - Configuration file to read, see below for existing list
  local port  - Optional local port number to use for the tunnel to avoid conflicts

EOM

  print_available_configs
}

[[ $# -gt 3 ]] && usage && exit 1

action="$1"
if [[ $action == "stop" ]]; then
    egrep '^LOCAL_TUNNEL_PORT=' ${CONFIG_DIR}/*.config | awk -F= '{print $2}' | sort -u | while read port; do
      echo "Stopping existing tunnel on local port $port"
      cleanup $port
    done
    echo "Tunnels stopped"
elif [[ $action == "start" ]]; then
  [[ $# -ne 2 && $# -ne 3 ]] && usage && exit 1
  config_file="${CONFIG_DIR}/${2}.config"
  [[ ! -f $config_file ]] && echo "Unknown configuration \"$1\": no file [$(basename ${config_file})] in [${CONFIG_DIR}]" && print_available_configs && exit 1

  source $config_file

  for var in PROJECT CLUSTER POD_NAME REMOTE_APPLICATION_PORT LOCAL_APPLICATION_PORT LOCAL_TUNNEL_PORT; do
    [[ -z ${!var} ]] && echo "$var must be specified in [$config_file]" && exit 1
  done

  if [[ "$PROJECT" == "actin-emc" ]]; then
    [[ -z ${NAMESPACE} ]] && echo "NAMESPACE must be specified in [$config_file] when PROJECT is actin-emc" && exit 1
  fi

  namespace=""
  [[ -n "$NAMESPACE" ]] && namespace="-n $NAMESPACE"

  GC="gcloud --project $PROJECT compute"
    echo "Fetching cluster credentials"
    gcloud container clusters get-credentials "$CLUSTER" --zone europe-west4 --project "$PROJECT"
    echo "Attempting to determine name of bastion VM"
    read instance zone <<< $($GC instances list | grep bastion | head -n1 | awk '{print $1, $2}')
    [[ -z $instance || -z $zone ]] && echo "Cannot locate bastion VM instance" && exit 1
    ps aux | grep gcloud | grep ssh | grep -- "-L $LOCAL_TUNNEL_PORT:localhost:$REMOTE_TUNNEL_PORT" >/dev/null
    if [[ $? -ne 0 ]]; then
      echo "Establishing tunnel to $instance on port $LOCAL_TUNNEL_PORT"
      $GC ssh $instance --tunnel-through-iap --zone $zone -- -L "$LOCAL_TUNNEL_PORT:localhost:$REMOTE_TUNNEL_PORT" -N -q -f
    else 
      echo "Re-using existing SSH tunnel to bastion"
    fi
    echo "Forwarding $REMOTE_APPLICATION_PORT to localhost:$LOCAL_APPLICATION_PORT for ${POD_NAME}"
    k8 $LOCAL_TUNNEL_PORT port-forward $(k8 $LOCAL_TUNNEL_PORT get pods $namespace | grep $POD_NAME | awk '{print $1}') $namespace \
        "$LOCAL_APPLICATION_PORT:$REMOTE_APPLICATION_PORT" &
    echo "Tunnel started; access $POD_NAME at http://localhost:$LOCAL_APPLICATION_PORT"
else
    echo "Unknown action \"$action\""
    usage
    exit 1
fi


