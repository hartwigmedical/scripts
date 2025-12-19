#!/usr/bin/env bash

REMOTE_TUNNEL_PORT=8080
KUBE_PORT=8888
CONFIG_DIR="$(readlink -f "$(dirname "$0")")/.tunnel_configurations"

function print_available_configs() {
    echo "Known configurations:"
    for c in "$(ls $CONFIG_DIR)"; do
        echo "${c}" | sed -e 's/^/  /' -e 's/\.config$//'
    done
}

function cleanup() {
    echo "Stopping all tunnels to bastion VMs.  This script does not support multiple tunnels running on one local machine."
    ps aux | grep "bastion" | grep 'start-iap-tunnel' | awk '{print $2}' | while read pid; do
        [[ -n $pid ]] && kill $pid 
        echo "  Process [$pid] killed"
    done
    echo "Cleanup complete"
}

function k8() {
    port=$1 && shift 
    HTTPS_PROXY=localhost:$port kubectl $@
}

function start_remote_proxy() {
    local gc_cli="$1"
    local instance="$2"
    local zone="$3"
    local namespace="$4"
    local project="$5"
    local cluster="$6"
    local services="$7"

    local remote_script
    remote_script=$(cat <<EOF
set -euo pipefail
LOG_FILE="\${HOME}/actin-ui-proxy.log"
PORT_LINE=""

get_port_line() {
    (grep "PROXY_PORT=" "\${LOG_FILE}" 2>/dev/null || true) | tail -n1
}

mkdir -p "\$(dirname "\${LOG_FILE}")"
[[ -f "\${LOG_FILE}" ]] || : > "\${LOG_FILE}"
PORT_LINE=\$(get_port_line)

if [[ -z "\${PORT_LINE}" ]]; then
    echo "Starting actin-ui-proxy for \${USER}"
    nohup actin-ui-proxy $project $namespace $cluster $services >> "\${LOG_FILE}" 2>&1 &
    for _ in \$(seq 1 30); do
        PORT_LINE=\$(get_port_line)
        if [[ -n "\${PORT_LINE}" ]]; then
            break
        fi
        sleep 1
    done
else
    echo "Found existing PROXY_PORT entry in \${LOG_FILE}; reusing running proxy"
fi

if [[ -z "\${PORT_LINE}" ]]; then
    echo "Failed to find PROXY_PORT in \${LOG_FILE}" >&2
    exit 1
fi

echo "\${PORT_LINE}"
EOF
)

    local proxy_output remote_port
    proxy_output=$($gc_cli ssh "$instance" --tunnel-through-iap --zone "$zone" --command "bash -lc '$remote_script'")
    remote_port=$(echo "$proxy_output" | awk -F= '/PROXY_PORT=/{print $2; exit}')
    if [[ -z "$remote_port" ]]; then
        echo "Unable to determine remote proxy port from actin-ui-proxy output" >&2
        echo "$proxy_output" >&2
        exit 1
    fi

    echo "$proxy_output" >&2
    echo "$remote_port"
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

cleanup

action="$1"
if [[ $action == "stop" ]]; then
    echo "Tunnels stopped"
elif [[ $action == "start" ]]; then
    [[ $# -ne 2 && $# -ne 3 ]] && usage && exit 1
    config_file="${CONFIG_DIR}/${2}.config"
    [[ ! -f $config_file ]] && echo "Unknown configuration \"$1\": no file [$(basename ${config_file})] in [${CONFIG_DIR}]" && print_available_configs && exit 1

    source $config_file

    for var in PROJECT BASTION PORT CONTEXT_ROOT CLUSTER NAMESPACE; do
        [[ -z ${!var} ]] && echo "$var must be specified in [$config_file]" && exit 1
    done

    GC="gcloud --project $PROJECT compute"
    echo "Attempting to determine name of bastion VM"
    read instance zone <<< $($GC instances list | grep $BASTION | head -n1 | awk '{print $1, $2}')
    [[ -z $instance || -z $zone ]] && echo "Cannot locate bastion VM instance" && exit 1
    ps aux | grep gcloud | grep ssh | grep -- "-L $PORT:localhost:$REMOTE_TUNNEL_PORT" >/dev/null
    if [[ $? -ne 0 ]]; then
        REMOTE_PROXY_PORT=$(start_remote_proxy "$GC" "$instance" "$zone" "$NAMESPACE" "$PROJECT" "$CLUSTER" "$CONTEXT_ROOT")
        echo "Using bastion proxy on port $REMOTE_PROXY_PORT"
        $GC ssh $instance --tunnel-through-iap --zone $zone -- -L "$REMOTE_PROXY_PORT:localhost:$REMOTE_TUNNEL_PORT" -L "$KUBE_PORT:localhost:$KUBE_PORT" -N -q -f
    else 
        echo "Re-using existing SSH tunnel to bastion"
    fi
    link="http://localhost:$PORT/$CONTEXT_ROOT"
    echo "Tunnel started; access $2 at: $link. Opening..."
    open "$link"
else
    echo "Unknown action \"$action\""
    usage
    exit 1
fi
