# Overview

Extract events from data on `datastore` for bulk loading by `turquoise`. Flow is:

1. An operator runs something that updates one of the files of interest
1. (time passes)
1. One of these scripts is run, checks for differences and generates a JSON file containing new events
1. The JSON files are pushed to a GCP bucket
1. A Turquoise instance is hopefully notified and collects the events from the bucket

To start from scratch and have all events re-published from all time, delete the relevant contents of 
the `data` directory and run again.

The `reported_date.sh` script is tied to the specific data format of the input file. Beyond that all processing
happens on Turquoise.

# Getting Turqoise to Pick Up Events

`gsutil notification create -t turquoise.storage -f json gs://your_bucket_name`

Make sure the Turquoise user has `Storage Object Viewer` and `Storage Legacy Bucket Reader` on the bucket.


