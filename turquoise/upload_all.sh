#!/usr/bin/env bash
#
# Wrapper to allow all scripts to be run from cron easily

/data/common/repos/scripts/turquoise/report_date.sh && \
/data/common/repos/scripts/turquoise/arrived_date.sh && \
gsutil cp ${HOME}/.turquoise/reported.json gs://turquoise-events-prod/ &&
gsutil cp ${HOME}/.turquoise/lims.json gs://turquoise-lims-prod/
