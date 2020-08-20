#!/usr/bin/env bash

#put bucket name (hmf-dr-XXX) at $1.

echo "[INFO] Files in the bucket" $1":"
gsutil du -h gs://$1/

echo "[INFO] Persmissions of the bucket"  $1":"
gsutil -u hmf-share iam get gs://$1/