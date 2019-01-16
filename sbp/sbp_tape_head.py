#!/bin/python
# check object tiering status - thoutenbos@schubergphilis.com

import re
import json
import sys
import requests
from dateutil import parser
from requests_aws4auth import AWS4Auth

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
  print 'Usage: head.py <bucket/object>'
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
    restore = parse_amz_restore(r.headers['x-amz-restore'])

    if restore['ongoing-request'] == 'true':
      print 'Object is being restored'
    else:
      print 'Object is restored till ' + str(restore['expiry-date'])
  
elif 'x-gmt-object-tieringinfo' in r.headers:
  print 'Object is tiered to Cold storage'

elif r.status_code != 200:
  print 'Object does not exist'

else:
  print 'Object is not tiered'
