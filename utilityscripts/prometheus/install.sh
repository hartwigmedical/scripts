####
install_base_dir=/opt/prometheus
data_base_dir=/data/prometheus
github_base_url=https://raw.githubusercontent.com/hartwigmedical/scripts/master/utilityscripts/prometheus
####

function install_prometheus() {
    add_user_and_group

    prometheus_version=1.4.1
    data_dir=${data_base_dir}/prometheus

    cd $(mktemp -d)
    prometheus_archive=prometheus-${prometheus_version}.linux-amd64.tar.gz
    wget -q https://github.com/prometheus/prometheus/releases/download/v${prometheus_version}/${prometheus_archive}

    prometheus_dir=${install_base_dir}/prometheus
    mkdir -p ${prometheus_dir}
    tar -zxf ${prometheus_archive} --strip-components=1 -C ${prometheus_dir}
    mkdir -p ${data_dir}
    chown -R prometheus.prometheus ${data_dir}

    service=prometheus.service
    cat << EOF > /usr/lib/systemd/system/${service}
[Unit]
Description=The Prometheus monitoring system and time series database.
Documentation=https://prometheus.io
After=network.target

[Service]
EnvironmentFile=-/etc/default/prometheus
User=prometheus
WorkingDirectory=${prometheus_dir}
ExecStart=${prometheus_dir}/prometheus \\
          -storage.local.path=${data_dir} \\
          -web.external-url=http://localhost:9090 \\
          \$PROMETHEUS_OPTS
ExecReload=/bin/kill -HUP \$MAINPID
Restart=on-failure

[Install]
WantedBy=multi-user.target
EOF
    systemd_enable ${service}

    wget -q ${github_base_url}/tunnel_node_exporters -O ${install_base_dir}/tunnel_node_exporters
    chmod +x ${install_base_dir}/tunnel_node_exporters

    cat <<EOF > ${prometheus_dir}/prometheus.yml
global:
  scrape_interval: 15s
  evaluation_interval: 15s
scrape_configs:
  - job_name: prometheus
    static_configs:
      - targets:
        - localhost:9090
  - job_name: node
    static_configs:
      - targets:
          - localhost:9100
        labels:
          instance: datastore
      - targets:
          - localhost:9101
        labels:
          instance: crunch001
      - targets:
          - localhost:9102
        labels:
          instance: crunch002
      - targets:
          - localhost:9103
        labels:
          instance: crunch003
rule_files:
  - node.rules
alerting:
  alertmanagers:
    - static_configs:
      - targets:
        - localhost:9093
EOF

    cat <<EOF > ${prometheus_dir}/node.rules
ALERT InstanceDown
  IF up == 0
  FOR 5m
  ANNOTATIONS {
    summary = "Instance {{ \$labels.instance }} down",
    description = "{{ \$labels.instance }} has been down for more than 5 minutes.",
  }
ALERT AlertingAlive
  IF day_of_month() <= 7 and day_of_week() == 1 and hour() == 10 and minute() >= 00 and minute() < 01
  ANNOTATIONS {
    summary = "Alerting is working",
    description = "Firing once per day, only to ensure alerting is working.",
  }
ALERT RAIDFailure
  IF megaraid_failed > 0
  ANNOTATIONS {
    summary = "{{ \$labels.instance }} RAID failure",
    description = "{{ \$labels.instance }} has a RAID failure.",
  }
ALERT RAIDDegraded
  IF megaraid_degraded > 0
  ANNOTATIONS {
    summary = "{{ \$labels.instance }} RAID degradation",
    description = "{{ \$labels.instance }} has a RAID degradation.",
  }
ALERT StorCLINotCompleted
  IF min(time() - storcli_completion_time{job="node"}) > 60 * 60
  FOR 10m
  ANNOTATIONS {
    summary = "{{ \$labels.instance }} storcli not running",
    description = "{{ \$labels.instance }} storcli has not completed successfully in over an hour.",
  }
ALERT DiskWillFillIn1Hour
  IF predict_linear(node_filesystem_free{job='node'}[1h], 2 * 60 * 60) < 0
  FOR 1h
  ANNOTATIONS {
    summary = "{{ \$labels.instance }} disk filling",
    description = "{{ \$labels.instance }} disk is predicted to fill within 1 hour.",
  }
EOF
}

function install_node_exporter() {
    add_user_and_group

    node_exporter_version=0.13.0
    data_dir=${data_base_dir}/node_exporter

    cd $(mktemp -d)
    node_exporter_archive=node_exporter-${node_exporter_version}.linux-amd64.tar.gz
    wget -q https://github.com/prometheus/node_exporter/releases/download/v${node_exporter_version}/${node_exporter_archive}
    node_exporter_dir=${install_base_dir}/node_exporter
    mkdir -p ${node_exporter_dir}
    tar -zxf ${node_exporter_archive} --strip-components=1 -C ${node_exporter_dir}
    mkdir -p ${data_dir}
    chown -R prometheus.prometheus ${data_dir}

    service=node_exporter.service
    cat << EOF > /usr/lib/systemd/system/${service}
[Unit]
Description=Prometheus exporter for machine metrics, written in Go with pluggable metric collectors.
Documentation=https://github.com/prometheus/node_exporter
After=network.target

[Service]
EnvironmentFile=-/etc/default/node_exporter
User=prometheus
WorkingDirectory=${node_exporter_dir}
ExecStart=${node_exporter_dir}/node_exporter \\
    --collector.textfile.directory=${data_dir} \\
    \$NODE_EXPORTER_OPTS
Restart=on-failure

[Install]
WantedBy=multi-user.target
EOF
    systemd_enable ${service}

    wget -q ${github_base_url}/storcli.py -O ${install_base_dir}/storcli.py
    chmod +x ${install_base_dir}/storcli.py
    cat << EOF > /etc/cron.d/prometheus-storcli
*/5 * * * * root ${install_base_dir}/storcli.py > ${data_dir}/storcli.prom.\$\$ && echo storcli_completion_time \$(date +\%s) >> ${data_dir}/storcli.prom.\$\$ && mv ${data_dir}/storcli.prom.\$\$ ${data_dir}/storcli.prom
EOF
}

function install_alertmanager() {
    if [ $# -lt 1 ];
    then
	echo "supply a webhook url from slack custom integrations"
	return
    fi
    web_hook_url=$1 && shift

    add_user_and_group

    alertmanager_version=0.5.1
    data_dir=${data_base_dir}/alertmanager

    cd $(mktemp -d)
    alertmanager_archive=alertmanager-${alertmanager_version}.linux-amd64.tar.gz
    wget -q https://github.com/prometheus/alertmanager/releases/download/v${alertmanager_version}/${alertmanager_archive}
    alertmanager_dir=${install_base_dir}/alertmanager
    mkdir -p ${alertmanager_dir}
    tar -zxf ${alertmanager_archive} --strip-components=1 -C ${alertmanager_dir}
    mkdir -p ${data_dir}
    chown -R prometheus.prometheus ${data_dir}

    service=alertmanager.service
    cat << EOF > /usr/lib/systemd/system/${service}
[Unit]
Description=Prometheus Alertmanager.
Documentation=https://github.com/prometheus/alertmanager
After=network.target

[Service]
EnvironmentFile=-/etc/default/alertmanager
User=prometheus
WorkingDirectory=${alertmanager_dir}
ExecStart=${alertmanager_dir}/alertmanager \\
          -storage.path=${data_dir} \\
          \$ALERTMANAGER_OPTS
ExecReload=/bin/kill -HUP \$MAINPID
Restart=on-failure

[Install]
WantedBy=multi-user.target
EOF
    systemd_enable ${service}

    cat <<EOF > ${alertmanager_dir}/alertmanager.yml
route:
 receiver: slack_general
 routes:
  - match:
      severity: slack
    receiver: slack_general

receivers:
- name: slack_general
  slack_configs:
  - api_url: ${web_hook_url}
    channel: '#alerts'
    send_resolved: true
EOF
}

function add_user_and_group() {
    getent group prometheus >/dev/null || sudo groupadd -r prometheus
    getent passwd prometheus >/dev/null || \
	useradd -r -g prometheus -d ${install_base_dir} -s /sbin/nologin \
            -c "Prometheus services" prometheus
}

function systemd_enable() {
    systemctl daemon-reload
    systemctl enable $1
}
