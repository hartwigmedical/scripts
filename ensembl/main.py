#!/usr/bin/env python3
import json
import logging
import re
import sys
import time
from concurrent.futures.thread import ThreadPoolExecutor
from copy import deepcopy
from threading import Lock
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

# set up logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.DEBUG,
    datefmt='%Y-%m-%d %H:%M:%S')


V37 = "37"
V38 = "38"

class StringCollector(object):
    def __init__(self):
        self.strings = []

    def get_all(self):
        return deepcopy(self.strings)

    def add(self, string):
        self.strings.append(string)


class BaseRestClient(object):
    def __init__(self, server, reqs_per_sec, max_wait_time=30):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.max_wait_time = max_wait_time
        self._req_count = 0
        self._last_req = 0
        self._wait_lock = Lock()
        self._wait_time_owed = 0

    def perform_rest_action(self, endpoint, headers=None, params=None):
        self._wait_lock.acquire()
        self._do_rate_limiting()

        self._wait_lock.release()

        if headers is None:
            headers = {}

        if 'Content-Type' not in headers:
            headers['Content-Type'] = 'application/json'

        if params:
            request_url = self.server + endpoint + '?' + urlencode(params)
        else:
            request_url = self.server + endpoint

        data = None

        try:
            request = Request(request_url, headers=headers)
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
                    logging.info("Retry after: {retry}".format(retry=retry))
                    self._wait_time_owed += float(retry) + 1
                    return self.perform_rest_action(endpoint, headers, params)
            elif e.code == 503:
                logging.error("Maybe recoverable error for request: " + request_url + "\nerror: " + repr(e))
                self._wait_time_owed += 5
                return self.perform_rest_action(endpoint, headers, params)
            else:
                logging.error("Fatal error for request: " + request_url + "\nerror: " + repr(e))
                raise ValueError(
                    ('Request failed for {request_url}: Status code: {e.code} Reason: {e.reason},\n'
                     'Full error: {full_error}').format(request_url=request_url, e=e, full_error=repr(e))
                )
        except URLError as e:
            logging.error("Maybe recoverable URL error for request: " + request_url + "\nerror: " + repr(e))
            return self.perform_rest_action(endpoint, headers, params)
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
            # enforce maximum wait time
            intended_wait_time = self._wait_time_owed
            actual_sleep_time = min(intended_wait_time, self.max_wait_time)
            logging.debug(
                "Sleep {actual_sleep_time} seconds, reduce owed wait by {intended_wait_time} seconds".format(
                    actual_sleep_time=actual_sleep_time, intended_wait_time=intended_wait_time
                )
            )
            time.sleep(actual_sleep_time)
            self._wait_time_owed -= intended_wait_time


class EnsemblRestClient(object):

    FULL_CIGAR_MATCH = re.compile(r"\d+M")

    def __init__(self, version: str, reqs_per_sec=15):
        if version == V38:
            self._rest_client = BaseRestClient('https://rest.ensembl.org', reqs_per_sec)
        elif version == V37:
            self._rest_client = BaseRestClient('https://grch37.rest.ensembl.org', reqs_per_sec)
        else:
            raise ValueError(f"Unrecognized ref genome version number: {version}")

    def get_coordinate_to_translated_position_and_sequences(
            self, species, coordinates, source_coordinate_system, target_coordinate_system, 
            warning_collector, max_workers=None
    ):
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            coordinate_to_future = {
                coordinate: executor.submit(
                    self.get_translated_position_and_sequences, species, coordinate,
                    source_coordinate_system, target_coordinate_system, warning_collector
                )
                for coordinate in coordinates
            }
            coordinate_to_result = {}
            for coordinate, future in coordinate_to_future.items():
                try:
                    coordinate_to_result[coordinate] = future.result()
                except Exception as e:
                    warning_collector.add(e)
            return coordinate_to_result

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
            external_names = {
                overview["external_name"] for overview in gene_overviews_for_id if "external_name" in overview.keys()
            }
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

    def get_translated_position_and_sequences(
            self, species, coordinate, source_coordinate_system, target_coordinate_system, warning_collector):
        logging.info("Start with getting translation for coordinate {coordinate}".format(coordinate=coordinate))
        sequence_pad_size = 5

        chrom, position = coordinate
        translated_position = self._request_translated_position(
            species, chrom, position, source_coordinate_system, target_coordinate_system, warning_collector
        )
        target_sequence = self._request_sequence(
            species, chrom, translated_position - sequence_pad_size, translated_position + sequence_pad_size,
            target_coordinate_system)
        source_sequence = self._request_sequence(
            species, chrom, position - sequence_pad_size, position + sequence_pad_size, source_coordinate_system)
        return translated_position, source_sequence, target_sequence

    def get_ensembl_id_to_nm_transcript_names(self, species, ensembl_ids, warning_collector, max_workers=None):
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_symbol = {
                ensembl_id: executor.submit(self.get_nm_transcript_names, species, ensembl_id, warning_collector)
                for ensembl_id in ensembl_ids
            }
        return {ensembl_id: future.result() for ensembl_id, future in future_to_symbol.items()}

    def get_nm_transcript_names(self, species, ensembl_id, warning_collector):
        logging.info("Start with getting overview for ensembl_id {ensembl_id}".format(ensembl_id=ensembl_id))
        return self._request_nm_transcript_ids(species, ensembl_id)

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

    def _request_translated_position(self, species, chrom, pos, source_coordinate_system, target_coordinate_system, warning_collector):
        try:
            try:
                target_start, target_end = self._request_translated_range(
                    species, chrom, pos, pos, source_coordinate_system, target_coordinate_system, warning_collector
                )

                assert target_start == target_end, "Translated start and end are different"
                return target_start
            except AssertionError as e:
                warning_collector.add(str(e))
                try:
                    offset = 10
                    adjusted_source_start = pos - offset
                    adjusted_source_end = pos + offset
                    adjusted_target_start, adjusted_target_end = self._request_translated_range(
                        species, chrom, adjusted_source_start, adjusted_source_end,
                        source_coordinate_system, target_coordinate_system, warning_collector
                    )
                    assert adjusted_target_start + offset == adjusted_target_end - offset, \
                        "Estimated translated start and end are different"
                    return adjusted_target_start + offset
                except AssertionError as e:
                    warning_collector.add(str(e))
                    try:
                        offset = 30
                        adjusted_source_start = pos - offset
                        adjusted_source_end = pos + offset
                        adjusted_target_start, adjusted_target_end = self._request_translated_range(
                            species, chrom, adjusted_source_start, adjusted_source_end,
                            source_coordinate_system, target_coordinate_system, warning_collector
                        )
                        assert adjusted_target_start + offset == adjusted_target_end - offset, \
                            "Second estimated translated start and end are different"
                        return adjusted_target_start + offset
                    except AssertionError as e:
                        warning_collector.add(str(e))
                        try:
                            offset = 100
                            adjusted_source_start = pos - offset
                            adjusted_source_end = pos + offset
                            adjusted_target_start, adjusted_target_end = self._request_translated_range(
                                species, chrom, adjusted_source_start, adjusted_source_end,
                                source_coordinate_system, target_coordinate_system, warning_collector
                            )
                            assert adjusted_target_start + offset == adjusted_target_end - offset, \
                                "Third estimated translated start and end are different"
                            return adjusted_target_start + offset
                        except AssertionError as e:
                            warning_collector.add(str(e))
                            try:
                                offset = 300
                                adjusted_source_start = pos - offset
                                adjusted_source_end = pos + offset
                                adjusted_target_start, adjusted_target_end = self._request_translated_range(
                                    species, chrom, adjusted_source_start, adjusted_source_end,
                                    source_coordinate_system, target_coordinate_system, warning_collector
                                )
                                assert adjusted_target_start + offset == adjusted_target_end - offset, \
                                    "Fourth estimated translated start and end are different"
                                return adjusted_target_start + offset
                            except AssertionError as e:
                                warning_collector.add(str(e))
        except Exception as e:
            warning_collector.add(str(e))

    def _request_translated_range(self, species, chrom, start, end, source_coordinate_system, target_coordinate_system,
                                  warning_collector):
        answer = self._rest_client.perform_rest_action(
            endpoint="/map/{0}/{1}/{2}:{3}..{4}:1/{5}".format(
                species, source_coordinate_system, chrom, start, end, target_coordinate_system),
            headers={"Content-Type": "application/json"}
        )
        mappings = answer["mappings"]
        assert mappings, "No mappings available: {chrom}:{start}-{end}\nanswer={answer}".format(
            chrom=chrom, start=start, end=end, answer=answer)
        assert len(mappings) == 1, "More than one mapping: {chrom}:{start}-{end}\nanswer={answer}".format(
            chrom=chrom, start=start, end=end, answer=answer)
        source_answer = mappings[0]["original"]
        target_answer = mappings[0]["mapped"]
        assert source_answer["assembly"] == source_coordinate_system, "Incorrect source coord system: {chrom}:{start}-{end}\nanswer={answer}".format(
            chrom=chrom, start=start, end=end, answer=answer)
        assert source_answer["seq_region_name"] == chrom, "Incorrect source chrom: {chrom}:{start}-{end}\nanswer={answer}".format(
            chrom=chrom, start=start, end=end, answer=answer)
        assert source_answer["strand"] == 1, "Incorrect source strand: {chrom}:{start}-{end}\nanswer={answer}".format(
            chrom=chrom, start=start, end=end, answer=answer)
        assert target_answer["assembly"] == target_coordinate_system, "Incorrect target coord system: {chrom}:{start}-{end}\nanswer={answer}".format(
            chrom=chrom, start=start, end=end, answer=answer)
        assert target_answer["seq_region_name"] == chrom, "Incorrect target chrom: {chrom}:{start}-{end}\nanswer={answer}".format(
            chrom=chrom, start=start, end=end, answer=answer)

        if target_answer["strand"] != 1:
            warning_collector.add("Incorrect target strand: {chrom}:{start}-{end}\nanswer={answer}".format(
            chrom=chrom, start=start, end=end, answer=answer))

        if source_answer["start"] == start and source_answer["end"] == end:
            return target_answer["start"], target_answer["end"]
        else:
            estimated_real_target_start = target_answer["start"] - (source_answer["start"] - start)
            estimated_real_target_end = target_answer["end"] + (end - source_answer["end"])
            warning = "\n".join([
                "Could not translate intended range precisely:",
                "Intended source range: {chrom}:{start}-{end}".format(chrom=chrom, start=start, end=end),
                "Realized source range: {source_seq_region_name}:{source_start}-{source_end}".format(
                    source_seq_region_name=source_answer["seq_region_name"],
                    source_start=source_answer["start"],
                    source_end=source_answer["end"],
                ),
                "Realized target range: {target_seq_region_name}:{target_start}-{target_end}".format(
                    target_seq_region_name=target_answer["seq_region_name"],
                    target_start=target_answer["start"],
                    target_end=target_answer["end"],
                ),
                "Estimated real target range: {target_seq_region_name}:{target_start}-{target_end}".format(
                    target_seq_region_name=target_answer["seq_region_name"],
                    target_start=estimated_real_target_start,
                    target_end=estimated_real_target_end,
                ),
                "Estimation potentially correct: {correct}".format(
                    correct=(estimated_real_target_end - estimated_real_target_start == end - start))
            ])

            warning_collector.add(warning)
            return estimated_real_target_start, estimated_real_target_end

    def _request_sequence(self, species, chrom, start, end, coordinate_system):
        answer = self._rest_client.perform_rest_action(
            endpoint="/sequence/region/{0}/{1}:{2}..{3}:1".format(
                species, chrom, start, end),
            headers={"Content-Type": "application/json"},
            params={"coord_system_version": coordinate_system}
        )
        return answer["seq"]

    def _request_ensembl_ids(self, species, symbol):
        answers = self._rest_client.perform_rest_action(
            endpoint="/xrefs/symbol/{0}/{1}".format(species, symbol),
            params={"object_type": "gene"}
        )
        return {answer["id"] for answer in answers if answer["id"].startswith("ENSG")}

    def _request_nm_transcript_ids(self, species, ensembl_id):
        if ensembl_id[:4] == "ENSG":
            object_type = "gene"
        elif ensembl_id[:4] == "ENST":
            object_type = "transcript"
        else:
            raise ValueError("Unrecognized Ensembl ID type")

        transcript_ids = list()
        crossrefs = self._rest_client.perform_rest_action(
            '/xrefs/id/{0}'.format(ensembl_id),
            params={"species": species, "object_type": object_type, "all_levels": 1}
        )
        for crossref in crossrefs:
            if crossref["db_display_name"] == "RefSeq mRNA" and ("cigar_line" not in crossref.keys() or self.FULL_CIGAR_MATCH.fullmatch(crossref["cigar_line"])):
                transcript_ids.append(crossref["primary_id"])
        return transcript_ids

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


def print_gene_id_and_name(species, symbols, output_file, version):
    client = EnsemblRestClient(version, reqs_per_sec=5)
    warning_collector = StringCollector()

    symbol_to_gene_overview = client.get_symbol_to_gene_overview(species, symbols, warning_collector, 15)

    symbol_overview_pairs = sorted([(symbol, overview) for symbol, overview in symbol_to_gene_overview.items()])
    with open(output_file, "w") as out_f:
        for symbol, overview in symbol_overview_pairs:
            if overview is None:
                out_f.write("\t".join([symbol, "NO_MATCH", symbol]) + "\n")
            else:
                out_f.write("\t".join([symbol, overview["id"], overview["external_name"]]) + "\n")

    warnings = warning_collector.get_all()
    if warnings:
        for warning in warnings:
            logging.error(warning)
        raise RuntimeError(f"Errors detected: {warnings}")


def print_nm_transcript_ids(species, gene_name_ensembl_id_tuples, output_file, version):
    client = EnsemblRestClient(version, reqs_per_sec=5)
    warning_collector = StringCollector()
    ensembl_ids = [pair[1] for pair in gene_name_ensembl_id_tuples]
    ensembl_id_to_nm_transcript_names = client.get_ensembl_id_to_nm_transcript_names(species, ensembl_ids, warning_collector, 15)

    with open(output_file, "w") as out_f:
        for gene_name, ensembl_id in gene_name_ensembl_id_tuples:
            out_f.write("\t".join([gene_name, ensembl_id, ";".join(ensembl_id_to_nm_transcript_names[ensembl_id])]) + "\n")

    warnings = warning_collector.get_all()
    if warnings:
        for warning in warnings:
            logging.error(warning)
        raise RuntimeError(f"Errors detected: {warnings}")


def determine_gene_ids_and_canonical_names(input_file, output_file, version):
    species = "human"
    with open(input_file) as f:
        symbols = f.read().splitlines()
    logging.debug((", ".join([symbol for symbol in symbols])).encode("ascii"))
    print_gene_id_and_name(species, symbols, output_file, version)


def get_coordinate_to_translated_position_and_sequences(
        species, coordinates, source_coordinate_system, target_coordinate_system, version):
    client = EnsemblRestClient(version, reqs_per_sec=5)
    warning_collector = StringCollector()
    range_to_translated_region = client.get_coordinate_to_translated_position_and_sequences(
        species, coordinates, source_coordinate_system, target_coordinate_system, warning_collector
    )
    return range_to_translated_region, warning_collector


def translate_coordinates(input_file, output_file, version):
    species = "human"
    with open(input_file) as f:
        lines = f.read().splitlines()

    source_coordinate_system, target_coordinate_system = lines[0].split("\t")
    coordinates = [tuple(line.split("\t")) for line in lines[1:]]
    cast_coordinates = [(chrom, int(pos)) for chrom, pos in coordinates]

    coordinate_to_translated_position_and_sequences, warning_collector = get_coordinate_to_translated_position_and_sequences(
        species, cast_coordinates, source_coordinate_system, target_coordinate_system, version
    )
    with open(output_file, "w") as out_f:
        for coord, translation in sorted(coordinate_to_translated_position_and_sequences.items(), key=lambda pair: pair[0]):
            chrom, pos = coord
            translated_pos, source_seq, target_seq = translation
            out_f.write("\t".join([chrom, str(pos), str(translated_pos), source_seq, target_seq]) + "\n")

    warnings = warning_collector.get_all()
    if warnings:
        for warning in warnings:
            logging.error(warning)
        raise RuntimeError(f"Errors detected: {warnings}")


def determine_nm_transcript_ids(input_file, output_file, version):
    species = "human"
    with open(input_file) as f:
        lines = f.read().splitlines()
    gene_name_ensembl_id_tuples = [line.split("\t") for line in lines]
    print_nm_transcript_ids(species, gene_name_ensembl_id_tuples, output_file, version)


def main(request_type, input_file, output_file, version):
    if request_type == "c":
        determine_gene_ids_and_canonical_names(input_file, output_file, version)
    elif request_type == "t":
        translate_coordinates(input_file, output_file, version)
    elif request_type == "n":
        determine_nm_transcript_ids(input_file, output_file, version)
    else:
        raise SyntaxError("Unknown type")


if __name__ == '__main__':
    if len(sys.argv) == 5:
        request_type_arg = sys.argv[1]
        input_file_arg = sys.argv[2]
        output_file_arg = sys.argv[3]
        version_arg = sys.argv[4]
    else:
        raise SyntaxError("Requires four arguments: type_of_query, input_file, output_file and reference genome version")
    if version_arg != V37 and version_arg != V38:
        raise SyntaxError(f"Ref genome version needs to be either {V37} or {V38}")

    main(request_type_arg, input_file_arg, output_file_arg, version_arg)
