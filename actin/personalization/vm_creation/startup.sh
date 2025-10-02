#!/bin/bash
set -e

# --- Install Docker ---
apt-get update
apt-get install -y apt-transport-https ca-certificates curl gnupg lsb-release 

curl -fsSL https://download.docker.com/linux/debian/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/debian \
  $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null

apt-get update
apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin

apt-get update
apt-get install -y pre-commit git


# --- Create developer user and add to docker group ---
useradd -m -s /bin/bash developer || true
usermod -aG docker developer || true

# --- Configure subuid/subgid ---
dev_uid=$(id -u developer)
dev_gid=$(id -g developer)
sed -i "/^developer:/d" /etc/subuid || true
sed -i "/^developer:/d" /etc/subgid || true
echo "developer:${dev_uid}:65536" >> /etc/subuid
echo "developer:${dev_gid}:65536" >> /etc/subgid

# --- Configure Docker daemon.json for userns-remap ---
mkdir -p /etc/docker
cat >/etc/docker/daemon.json <<EOF
{
  "userns-remap": "developer"
}
EOF

systemctl restart docker

mkdir /data
chown developer:developer /data
chmod 755 /data

rm /etc/motd
