#!/usr/bin/python
# Copyright 2018 Google, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Script is a slightly modified version of:
# https://github.com/GoogleCloudPlatform/python-docs-samples/blob/master/storage/signed_urls/generate_signed_urls.py

import argparse
import binascii
import collections
import datetime
import hashlib
import sys
import six
from six.moves.urllib.parse import quote
from google.oauth2 import service_account

def generate_url(account_file, bucket_object, expiration):
    http_method='GET'
    query_parameters=None
    headers=dict()
    headers['host'] = 'storage.googleapis.com'

    if expiration > 604800:
        print('Expiration Time can\'t be longer than 604800 seconds (7 days).')
        sys.exit(1)

    object_path = quote(six.ensure_binary(bucket_object), safe=b'/~')
    canonical_uri = '/{}'.format(object_path)

    datetime_now = datetime.datetime.utcnow()
    request_timestamp = datetime_now.strftime('%Y%m%dT%H%M%SZ')
    datestamp = datetime_now.strftime('%Y%m%d')

    google_credentials = service_account.Credentials.from_service_account_file(account_file)
    client_email = google_credentials.service_account_email
    credential_scope = '{}/auto/storage/goog4_request'.format(datestamp)
    credential = '{}/{}'.format(client_email, credential_scope)

    canonical_headers = ''
    ordered_headers = collections.OrderedDict(sorted(headers.items()))
    for k, v in ordered_headers.items():
        lower_k = str(k).lower()
        strip_v = str(v).lower()
        canonical_headers += '{}:{}\n'.format(lower_k, strip_v)

    signed_headers = ''
    for k, _ in ordered_headers.items():
        lower_k = str(k).lower()
        signed_headers += '{};'.format(lower_k)
    signed_headers = signed_headers[:-1]  # remove trailing ';'

    if query_parameters is None:
        query_parameters = collections.OrderedDict()
    query_parameters['userProject'] = 'hmf-database'
    query_parameters['X-Goog-Algorithm'] = 'GOOG4-RSA-SHA256'
    query_parameters['X-Goog-Credential'] = credential
    query_parameters['X-Goog-Date'] = request_timestamp
    query_parameters['X-Goog-Expires'] = expiration
    query_parameters['X-Goog-SignedHeaders'] = signed_headers

    canonical_query_string = ''
    for k, v in query_parameters.items():
        encoded_k = quote(str(k), safe='')
        encoded_v = quote(str(v), safe='')
        canonical_query_string += '{}={}&'.format(encoded_k, encoded_v)
    canonical_query_string = canonical_query_string[:-1]  # remove trailing ';'
    canonical_request = '\n'.join([http_method, canonical_uri, canonical_query_string, canonical_headers, signed_headers, 'UNSIGNED-PAYLOAD'])
    canonical_request_hash = hashlib.sha256(canonical_request.encode()).hexdigest()
    string_to_sign = '\n'.join(['GOOG4-RSA-SHA256', request_timestamp, credential_scope, canonical_request_hash])
    signature = binascii.hexlify(google_credentials.signer.sign(string_to_sign)).decode()
    host_name = 'https://storage.googleapis.com'
    signed_url = '{}{}?{}&X-Goog-Signature={}'.format(host_name, canonical_uri, canonical_query_string, signature)
    return signed_url

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('account_file', help='Path to your Google service account.')
    parser.add_argument('bucket_object', help='Your Cloud Storage bucket/object path.')
    parser.add_argument('expiration', help='Expiration Time.')

    args = parser.parse_args()
    url = generate_url(account_file=args.account_file, bucket_object=args.bucket_object, expiration=int(args.expiration))

    print(url)
