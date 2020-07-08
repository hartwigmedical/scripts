# Overview 

Scripts to extract events from data on `datastore` for bulk loading by `turquoise`. Flow is:

1. An operator runs something that updates one of the files of interest
1. (time passes)
1. One of these scripts is run, checks for differences and generates a JSON file containing new events
1. The JSON files are pushed to a GCP bucket
1. A Turquoise instance is hopefully notified and collects the events from the bucket

To start from scratch and have all events re-published from all time, delete the relevant contents of 
the `data` directory and run again.

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

# Getting Turqoise to Pick Up Events

`gsutil notification create -t turquoise.storage -f json gs://your_bucket_name`


