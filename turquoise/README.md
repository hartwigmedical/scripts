# Overview 

Script to extract events from LIMS data on `datastore` for bulk loading by `turquoise`. Flow is:

1. An operator runs an update script that results in the `/data/ops/lims/prod/lims.json` file being written
1. (time passes)
1. This script is run. It checks for differences between the current data and what was there last time and generates a JSON file
   for each event type containing all the new events.
1. The JSON files are pushed to a GCP bucket.

That's it for the local flow, then on GCP a Kubernetes `cron` job will pick up the JSON files and push their contents into the
BigQuery table.

As the script does its own `diff` it keeps track of all old data files. To start from scratch and have it re-publish all events
from the start of time, just delete the contents of the `data` directory wherever its working directory is and it will push
everything over again.

# Input Data

This script is highly tied to the input data format. The data needs to be well-formed JSON containing a `samples` object:

```
# jq '.samples' lims.json
...
  "FR03631783": {
    ...
    "arrival_date": "2016-04-18",
    ...
    "label": "CPCT",
    ...
    "report_date": "",
    "sample_name": "CPCT02020292T",
    "sample_source": "RNA-TISSUE",
    ...
  },
  ...
```

