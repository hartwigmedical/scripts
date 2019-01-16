#!/bin/python
# Request object restores - thoutenbos

import re
import sys
import json
import requests
from requests_aws4auth import AWS4Auth

if len(sys.argv) != 2:
  print 'Usage: restore.py <bucket/object>'
  sys.exit(0)

## SVL: get credentials on datastore
awsProfile = "download"
credFilePath = '/home/sbp/.aws/credentials'
with open(credFilePath) as fp:
    for line in fp:
        if re.match("\[" + awsProfile + "\]", line) is not None:
             accessKeyLine = fp.next().rstrip()
             secretKeyLine = fp.next().rstrip()
             accessKey = accessKeyLine.split(' = ')[1]
             secretKey = secretKeyLine.split(' = ')[1]
fp.close()

r = requests.post(
  'https://s3.object02.schubergphilis.com/' + sys.argv[1] + '?restore',
        auth=AWS4Auth(
            accessKey,
            secretKey,
            'eu-nl-prod01',
            's3'
       	),
       	timeout=30,
        data='<RestoreRequest><Days>7</Days></RestoreRequest>'
    )

print r.status_code
print r.headers
print r.text


