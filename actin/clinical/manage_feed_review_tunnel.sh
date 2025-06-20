#!/usr/bin/env bash

REMOTE_PORT=8888
REMOTE_APPLICATION_PORT=8080
APPLICATION_PORT=8082

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <PROJECT> <ACTION> [<PORT>]"
    echo "  <PROJECT>   - GCP project ID, e.g. 'actin-nki'"
    echo "  <ACTION>    - 'start' to start the tunnel, 'stop' to stop it"
    echo "  <PORT>      - Optional port number to use for the tunnel (default: $REMOTE_PORT)"
    echo "After starting the tunnel, you can access the feed review tool at http://localhost:$APPLICATION_PORT. HTTPS is not supported."
    exit 1
fi

PROJECT=$1
GC="gcloud --project $PROJECT compute"

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

action="$2"
[[ $action != "start" && $action != "stop" ]] && echo "Provide a verb (either \"start\" or \"stop\")" && exit 1

port="${3:-$REMOTE_PORT}"

if [[ $action == "stop" ]]; then
    echo "Stopping tunnel running on port $port"
    cleanup $port
    echo "Tunnel stopped"
else
    echo "Starting tunnel on port $port"
    cleanup $port
    read instance zone <<< $($GC instances list | grep bastion | head -n1 | awk '{print $1, $2}')
    [[ -z $instance || -z $zone ]] && echo "Cannot locate bastion VM instance" && exit 1
    $GC ssh $instance --tunnel-through-iap --zone $zone -- -L $port:localhost:$REMOTE_PORT -N -q -f
    echo "Forwarding $REMOTE_APPLICATION_PORT to localhost:$APPLICATION_PORT for feed review tool"
    k8 $port port-forward $(k8 $port get pods | grep actin-feed-review | awk '{print $1}') $APPLICATION_PORT:$REMOTE_APPLICATION_PORT &
    echo "Tunnel started - access the feed review tool at http://localhost:$APPLICATION_PORT"
fi
