#!/bin/python
# Request object restores - thoutenbos

import re
import sys
import json
import requests
import hashlib
import base64
from requests_aws4auth import AWS4Auth
from dateutil import parser
from datetime import timedelta, datetime
from pytz import timezone


def parse_amz_restore(input):
    # x-amz-restore: ongoing-request="true", expiry-date="Fri, 27 Apr 2018 08:49:28 GMT"
    result = dict()

    regex = re.compile(r'([^=]+)="([^"]+)"(, )?')
    for match in regex.finditer(input):
        key = match.group(1)
        val = match.group(2)

        if key == 'expiry-date':
            val = parser.parse(val)

        result[key] = val
    return result

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

r = requests.head(
  'https://s3.object02.schubergphilis.com/' + sys.argv[1],
        auth=AWS4Auth(
            accessKey,
            secretKey,
            'eu-nl-prod01',
            's3'
        ),
        timeout=30
    )

if 'x-amz-restore' in r.headers:
    parsed = parse_amz_restore(r.headers['x-amz-restore'])

    if 'expiry-date' in parsed and parsed['expiry-date'] > datetime.now(timezone('Europe/Amsterdam')) + timedelta(days=7):
        print 'Object already restored till ' + str(parsed['expiry-date']) + ' so not requesting new restore'
        sys.exit(0)

data='<RestoreRequest><Days>7</Days></RestoreRequest>'

r = requests.post(
  'https://s3.object02.schubergphilis.com/' + sys.argv[1] + '?restore',
        auth=AWS4Auth(
            accessKey,
            secretKey,
            'eu-nl-prod01',
            's3'
       	),
        headers={
            'Content-MD5': base64.b64encode(hashlib.md5(data).digest())
        },
       	timeout=30,
        data=data
    )

print r.status_code
print r.headers
print r.text


