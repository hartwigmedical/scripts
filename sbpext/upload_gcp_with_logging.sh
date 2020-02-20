#!/bin/bash

/data/repos/scripts/sbpext/upload_gcp.sh $1 $2 > $2/logs/upload_$(date +%Y%m%d%H%M) 2>&1