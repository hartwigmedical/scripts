#!/usr/bin/env python
import json
import logging
import sys
import time
from concurrent.futures.thread import ThreadPoolExecutor
from copy import deepcopy
from threading import Lock
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

# set up logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.DEBUG,
    datefmt='%Y-%m-%d %H:%M:%S')


class StringCollector(object):
    def __init__(self):
        self.strings = []

    def get_all(self):
        return deepcopy(self.strings)

    def add(self, string):
        self.strings.append(string)


class BaseRestClient(object):
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self._req_count = 0
        self._last_req = 0
        self._wait_lock = Lock()
        self._wait_time_owed = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        self._wait_lock.acquire()
        self._do_rate_limiting()

        self._wait_lock.release()

        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            request_url = self.server + endpoint + '?' + urlencode(params)
        else:
            request_url = self.server + endpoint

        data = None

        try:
            request = Request(request_url, headers=hdrs)
            logging.debug("Send request to url: {request_url}".format(request_url=request_url))
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                logging.error("Recoverable error for request: " + request_url + "\nerror: " + repr(e))
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    self._wait_time_owed += float(retry) + 1
                    return self.perform_rest_action(endpoint, hdrs, params)
            elif e.code == 503:
                logging.error("Maybe recoverable error for request: " + request_url + "\nerror: " + repr(e))
                self._wait_time_owed += 5
                return self.perform_rest_action(endpoint, hdrs, params)
            else:
                logging.error("Fatal error for request: " + request_url + "\nerror: " + repr(e))
                raise ValueError(
                    ('Request failed for {request_url}: Status code: {e.code} Reason: {e.reason},\n'
                     'Full error: {full_error}').format(request_url=request_url, e=e, full_error=repr(e))
                )

        return data

    def _do_rate_limiting(self):
        # check if we need to rate limit ourselves
        if self._req_count >= self.reqs_per_sec:
            delta = time.time() - self._last_req
            if delta <= 1:
                self._wait_time_owed += 1
            self._last_req = time.time()
            self._req_count = 0
        self._req_count += 1

        # do waiting
        while self._wait_time_owed > 0:
            sleep_time = self._wait_time_owed
            logging.debug("Sleep {sleep_time} seconds".format(sleep_time=sleep_time))
            time.sleep(sleep_time)
            self._wait_time_owed -= sleep_time


class EnsemblRestClient(object):
    def __init__(self, reqs_per_sec=15):
        self._rest_client = BaseRestClient('http://rest.ensembl.org', reqs_per_sec)

    def get_symbol_to_gene_overview(self, species, symbols, warning_collector, max_workers=None):
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_symbol = {
                symbol: executor.submit(self.get_gene_overview, species, symbol, warning_collector)
                for symbol in symbols
            }
            return {symbol: future.result() for symbol, future in future_to_symbol.items()}

    def get_gene_overview(self, species, symbol, warning_collector):
        logging.info("Start with getting overview for symbol {symbol}".format(symbol=symbol))
        suggested_ensembl_ids = self._request_ensembl_ids(species, symbol)

        possible_gene_overviews = []
        for ensembl_id in suggested_ensembl_ids:
            gene_overviews_for_id = self._request_gene_overviews(species, ensembl_id)
            external_names = {overview["external_name"] for overview in gene_overviews_for_id}
            synonyms = self._request_synonyms(species, ensembl_id)

            if symbol in external_names or symbol in synonyms:
                possible_gene_overviews.extend(gene_overviews_for_id)

        gene_overview = self._select_unique_gene_overview(possible_gene_overviews, species, symbol, warning_collector)

        logging.info("Got gene overviews for symbol {symbol}".format(symbol=symbol))

        return gene_overview

    def get_variants(self, species, symbol, warning_collector):
        overview = self.get_gene_overview(species, symbol, warning_collector)
        if overview is None:
            return None
        else:
            return self._request_variants(species, overview["id"])

    def _select_unique_gene_overview(self, gene_overviews, species, symbol, warning_collector):
        if not gene_overviews:
            return None
        elif len(gene_overviews) == 1:
            return gene_overviews.pop()

        # so more than 1 relevant gene_overview

        precise_name_match_overviews = [
            overview for overview in gene_overviews if overview["external_name"] == symbol
        ]

        if len(precise_name_match_overviews) == 1:
            possibilities_string = ", ".join(
                [overview["id"] + ": " + overview["external_name"] for overview in gene_overviews]
            )
            choice_warning = (
                "Warning: selected overview with exact name match out of multiple possible choices:\n "
                "species={species}, symbol={symbol}, possibilities={possibilities}".format(
                    species=species, symbol=symbol, possibilities=possibilities_string)
            )
            warning_collector.add(choice_warning)
            return precise_name_match_overviews[0]

        # so no unique best choice
        gene_overviews_string = "\n".join([self._pretty_overview_string(overview) for overview in gene_overviews])

        if precise_name_match_overviews:
            # so multiple matches by exact name
            choice_warning = (
                "Error: Multiple gene overviews, multiple matches by exact name. "
                "species: {species}, symbol: {symbol}, gene_overviews_string:\n"
                "{gene_overviews_string}"
            ).format(species=species, symbol=symbol, gene_overviews_string=gene_overviews_string)
            warning_collector.add(choice_warning)
            return None
        else:
            # so no matches by exact name
            choice_warning = (
                "Error: Multiple gene overviews, no matches by exact name. "
                "species: {species}, symbol: {symbol}, gene_overviews_string:\n"
                "{gene_overviews_string}"
            ).format(species=species, symbol=symbol, gene_overviews_string=gene_overviews_string)
            warning_collector.add(choice_warning)
            return None

    def _request_ensembl_ids(self, species, symbol):
        answers = self._rest_client.perform_rest_action(
            endpoint="/xrefs/symbol/{0}/{1}".format(species, symbol),
            params={"object_type": "gene"}
        )
        return {answer["id"] for answer in answers if answer["id"].startswith("ENSG")}

    def _request_synonyms(self, species, ensembl_id):
        synonyms = set()
        crossrefs = self._rest_client.perform_rest_action(
            '/xrefs/id/{0}'.format(ensembl_id),
            params={"species": species, "object_type": "gene"}
        )
        for crossref in crossrefs:
            synonyms.update(crossref["synonyms"])
        return synonyms

    def _request_gene_overviews(self, species, ensembl_id):
        answers = self._rest_client.perform_rest_action(
            '/overlap/id/{0}'.format(ensembl_id),
            params={"feature": "gene", "species": species}
        )
        return [overview for overview in answers if self._overview_is_relevant(overview, ensembl_id)]

    def _request_variants(self, species, ensembl_id):
        return self._rest_client.perform_rest_action(
            '/overlap/id/{0}'.format(ensembl_id),
            params={"feature": "variation", "species": species}
        )

    @staticmethod
    def _overview_is_relevant(overview, ensembl_id):
        return overview["id"] == ensembl_id and not overview["seq_region_name"].startswith("CHR")

    @staticmethod
    def _pretty_overview_string(overview):
        key_to_value_strings = ["{key}: {value}".format(key=key, value=value) for key, value in overview.items()]
        return "{" + ", ".join(sorted(key_to_value_strings)) + "}"


def print_gene_id_and_name(species, symbols):
    client = EnsemblRestClient(reqs_per_sec=5)
    warning_collector = StringCollector()

    symbol_to_gene_overview = client.get_symbol_to_gene_overview(species, symbols, warning_collector, 15)

    symbol_overview_pairs = sorted([(symbol, overview) for symbol, overview in symbol_to_gene_overview.items()])
    for symbol, overview in symbol_overview_pairs:
        if overview is None:
            print("\t".join([symbol, "NO_MATCH", symbol]))
        else:
            print("\t".join([symbol, overview["id"], overview["external_name"]]))

    for warning in warning_collector.get_all():
        print(warning)


def print_variants(species, symbols):
    # This is pretty slow. Create an equivalent of EnsemblRestClient.get_symbol_to_gene_overview to make it faster.
    client = EnsemblRestClient(reqs_per_sec=5)
    warning_collector = StringCollector()
    for symbol in symbols:
        variants = client.get_variants(species, symbol, warning_collector)
        if variants:
            for v in variants:
                print("{seq_region_name}:{start}-{end}:{strand} ==> {id} ({consequence_type})".format(**v))

    for warning in warning_collector.get_all():
        print(warning)


def main():
    species = "human"

    if len(sys.argv) == 2:
        input_file = sys.argv[1]
    else:
        raise SyntaxError("Requires one argument: input_file")

    with open(input_file) as f:
        symbols = f.read().splitlines()

    logging.debug((", ".join([symbol for symbol in symbols])).encode("ascii"))
    # symbols = ['ABCB1']
    # symbols = ['BRAF']
    # symbols = ['ADGRB3']
    # symbols = ['BAI3']
    # symbols = ['AKT3']
    # symbols = ['ALB']
    # symbols = ['BABAM1']
    # symbols = ['H3F3B']
    # symbols = ['MLL2']
    #
    print_gene_id_and_name(species, symbols)


if __name__ == '__main__':
    main()
