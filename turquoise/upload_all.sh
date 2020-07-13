#!/usr/bin/env bash
#
# Wrapper to allow all scripts to be run from cron easily

/data/common/repos/scripts/turquoise/report_date.sh && \
gsutil cp ${HOME}/.turquoise/*.json gs://turquoise-events-prod/ && \
gsutil cp /data/ops/lims/prod/lims.json gs://turquoise-lims-prod/
