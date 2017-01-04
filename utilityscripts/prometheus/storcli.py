#!/usr/bin/env python

# Script to parse StorCLI's JSON output and expose
# MegaRAID health as Prometheus metrics.
#
# Tested against StorCLI 'Ver 1.14.12 Nov 25, 2014'.
#
# StorCLI reference manual:
# http://docs.avagotech.com/docs/12352476
#
# Advanced Software Options (ASO) not exposed as metrics currently.
#
# JSON key abbreviations used by StorCLI are documented in the standard command
# output, i.e.  when you omit the trailing 'J' from the command.

import json
import subprocess

METRIC_PREFIX = 'megaraid_'
METRIC_CONTROLLER_LABELS = '{{controller="{}", model="{}"}}'


def get_store_cli_json():
    storcli_cmd = ['/opt/MegaRAID/storcli/storcli64', 'show', 'all', 'J']
    proc = subprocess.Popen(storcli_cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    return proc.communicate()[0]


def main():
    data = json.loads(get_store_cli_json())

    # It appears that the data we need will always be present in the first
    # item in the Controllers array
    status = data['Controllers'][0]

    metrics = {
            'status_code': status['Command Status']['Status Code'],
            'controllers': status['Response Data']['Number of Controllers'],
            }

    for name, value in metrics.iteritems():
        print("{}{} {}".format(METRIC_PREFIX, name, value))

    for controller in status['Response Data']['System Overview']:
        controller_index = controller['Ctl']
        model = controller['Model']
        labels = METRIC_CONTROLLER_LABELS.format(controller_index, model)

        controller_metrics = {
                # FIXME: Parse dimmer switch options
                # 'dimmer_switch':          controller['DS'],

                'battery_backup_healthy':   int(controller['BBU'] == 'Opt'),
                'degraded':                 int(controller['Hlth'] == 'Dgd'),
                'drive_groups':             controller['DGs'],
                'emergency_hot_spare':      int(controller['EHS'] == 'Y'),
                'failed':                   int(controller['Hlth'] == 'Fld'),
                'healthy':                  int(controller['Hlth'] == 'Opt'),
                'physical_drives':          controller['PDs'],
                'ports':                    controller['Ports'],
                'scheduled_patrol_read':    int(controller['sPR'] == 'On'),
                'virtual_drives':           controller['VDs'],

                # Reverse StorCLI's logic to make metrics consistent
                'drive_groups_optimal':     int(controller['DNOpt'] == 0),
                'virtual_drives_optimal':   int(controller['VNOpt'] == 0),
                }

        for name, value in controller_metrics.iteritems():
            print("{}{}{} {}".format(METRIC_PREFIX, name, labels, value))

main()
