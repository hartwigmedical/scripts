#!/bin/bash

set -a  
source .env
set +a

VM_NAME=actin-developer-vm-$(whoami)
GIT_REPO="git@github.com:hartwigmedical/actin-personalization.git"

# verify all variables are set
: "${SSH_KEY_FILE:?Need to set SSH_KEY_FILE non-empty}"
: "${GIT_NAME:?Need to set GIT_NAME non-empty}"
: "${GIT_EMAIL:?Need to set GIT_EMAIL non-empty}"

# create ssh config entry
echo "
Host $VM_NAME
    User developer
    ForwardAgent yes
    IdentityFile $SSH_KEY_FILE
    IdentitiesOnly yes
    HostName $VM_NAME
    Port 1215
    ProxyCommand gcloud compute start-iap-tunnel  --zone europe-west4-a --project actin-research --listen-on-stdin %h 22
" > ~/.ssh/personalization_vm_config

# ensure ssh config file exists
if [ ! -f ~/.ssh/config ]; then
    echo "Creating ~/.ssh/config file"
    touch ~/.ssh/config
    chmod 600 ~/.ssh/config
fi

# remove old entry from known_hosts
sed -i.bak "/$VM_NAME/d" ~/.ssh/known_hosts

# check if entry is in main ssh and if not include in main ssh config
if ! grep -q "Include ~/.ssh/personalization_vm_config" ~/.ssh/config; then
    echo "Include ~/.ssh/personalization_vm_config" >> ~/.ssh/config
fi

while ! gcloud compute ssh developer@$VM_NAME --zone=europe-west4-a --tunnel-through-iap --ssh-key-file=$SSH_KEY_FILE --command "exit 0" --quiet 2>/dev/null; do
    echo "Waiting for VM to become available..."
    sleep 5
done
echo "VM is available!"

ssh -t -A $VM_NAME "bash -l -c 'gcloud auth login'"

ssh -t -A $VM_NAME << EOF
    git config --global user.name "$GIT_NAME"
    git config --global user.email "$GIT_EMAIL"
    ssh-keyscan github.com >> ~/.ssh/known_hosts

    if [ ! -d ~/actin-personalization ]; then
        git clone $GIT_REPO ~/actin-personalization/
    fi
    
    cd ~/actin-personalization
    pre-commit install

    if [ ! -d /data/actin-resources-private ]; then
        echo "Setting up private resources..."
        gcloud source repos clone actin-resources-private /data/actin-resources-private --project=actin-build
    fi

    echo "Setting up reference data..."
    mkdir -p /data/actin-personalization-reference-patient-data/
    gsutil -m rsync -r gs://actin-personalization-reference-patient-data /data/actin-personalization-reference-patient-data/
    mkdir -p /data/hmf-crunch-actin-ncr-dataset/
    gsutil -m rsync -r gs://hmf-crunch-actin-ncr-dataset /data/hmf-crunch-actin-ncr-dataset/
EOF

echo "Setup complete! You can now connect to your VM using: ssh $VM_NAME"
