#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import sys
import traceback
import zlib
import os
import math
import contextlib
import struct
import argparse


READ_SIZE = 4096
DECOMPRESSOR_OPTIONS = zlib.MAX_WBITS | 16  # max window, gzip header/trailer
GZIP_MAGIC = "\037\213"


def optimistic_decompress(chunks, num_chunks, no_cache):
    decompressor, bam_chunks, unused = init_decompressor()
    for chunk_num, chunk in chunks:
        try:
            bam_chunk = decompressor.decompress(unused + chunk)
            if bam_chunk:
                if no_cache:
                    yield bam_chunk
                else:
                    bam_chunks.append(bam_chunk)
        except Exception:
            print_exception(chunk_num, decompressor)
            decompressor, bam_chunks, unused = init_decompressor(chunks, decompressor)
            continue

        if decompressor.unconsumed_tail:
            print("unexpected unconsumed tail of length {}".format(len(decompressor.unconsumed_tail)), file=sys.stderr)

        unused = decompressor.unused_data
        if unused:
            # new member, yield and start fresh, keeping unused portion
            for bam_chunk in bam_chunks:
                yield bam_chunk
            decompressor, bam_chunks, _ = init_decompressor()

        if chunk_num % 1000 == 0:
            print_progress(chunk_num, num_chunks)

    # gzip module only does this at the end, not on new members
    bam_chunk = decompressor.flush()
    if bam_chunk:
        bam_chunks.append(bam_chunk)

    for bam_chunk in bam_chunks:
        yield bam_chunk

    print_progress(num_chunks, num_chunks)


def init_decompressor(chunks=None, decompressor=None):
    if chunks:
        unused = find_next_member(chunks, decompressor)
    else:
        unused = b""
    return zlib.decompressobj(DECOMPRESSOR_OPTIONS), [], unused


def print_exception(chunk_num, decompressor):
    error_offset = chunk_num * READ_SIZE - len(decompressor.unconsumed_tail)
    print(file=sys.stderr)
    print("decompression error at offset {} (unused data: {})".format(error_offset, len(decompressor.unused_data)), file=sys.stderr)
    traceback.print_exc(file=sys.stderr)


def find_next_member(chunks, decompressor):
    unused = b""
    next_part = decompressor.unconsumed_tail.find(GZIP_MAGIC)
    if next_part != -1:
        unused = decompressor.unconsumed_tail[next_part:]
    else:
        for chunk_num, chunk in chunks:
            next_part = chunk.find(GZIP_MAGIC)
            if next_part != -1:
                unused = chunk[next_part:]
                break
    return unused


def print_progress(chunk_num, num_chunks):
    if num_chunks is not None:
        percent_complete = chunk_num / num_chunks if num_chunks != 0 else 0
        print("Processed chunk {:{width}}/{:{width}} ({:.1%})".format(
            chunk_num,
            num_chunks,
            percent_complete,
            width=len(str(num_chunks))),
              end="\r",
              file=sys.stderr)
        if chunk_num == num_chunks:
            print(file=sys.stderr)
            print("Done.", file=sys.stderr)
        sys.stderr.flush()


def parse_bam(chunks, num_chunks, print_reads, all_reads):
    (magic, l_text), unused = unpacker("<4si", b"", chunks, num_chunks)
    (text, n_ref), unused = unpacker("<{}si".format(l_text), unused, chunks, num_chunks)
    for i in range(n_ref):
        (l_name,), unused = unpacker("<i", unused, chunks, num_chunks)
        (name, l_ref), unused = unpacker("<{}si".format(l_name), unused, chunks, num_chunks)
    try:
        prev_read_name = ""
        while True:
            fmt_to_read_name = "<3i2I4i"
            (block_size, refID, pos, bin_mq_nl, flag_nc, l_seq, next_refID, next_pos, tlen), unused = unpacker(fmt_to_read_name, unused, chunks, num_chunks)
            midparse = True
            # bin = bin_mq_nl >> 16
            # mapq = (bin_mq_nl & 0xff00) >> 8
            l_read_name = bin_mq_nl & 0xff
            # flag = flag_nc >> 16
            n_cigar_op = flag_nc & 0xffff
            fmt_read_name = "<{}s".format(l_read_name - 1)
            fmt_cigar = "<{}I".format(n_cigar_op)
            fmt_seq = "<{}B".format((l_seq + 1) // 2 + 1)
            fmt_qual = "<{}c".format(l_seq)
            (read_name, ), unused = unpacker(fmt_read_name, unused, chunks, num_chunks)
            # before the other parsing in case the rest goes wrong
            process_read(prev_read_name, read_name, refID, pos, print_reads, all_reads)
            cigar, unused = unpacker(fmt_cigar, unused, chunks, num_chunks)
            seq, unused = unpacker(fmt_seq, unused, chunks, num_chunks)
            qual, unused = unpacker(fmt_qual, unused, chunks, num_chunks)
            l_tags = (block_size
                      - struct.calcsize(fmt_to_read_name)
                      - struct.calcsize(fmt_read_name)
                      - struct.calcsize(fmt_cigar)
                      - struct.calcsize(fmt_seq)
                      - struct.calcsize(fmt_qual) + 4)
            tags, unused = unpacker("<{}c".format(l_tags), unused, chunks, num_chunks)
            midparse = False
            prev_read_name = read_name
    except StopIteration:
        if midparse:
            raise
        else:
            print_progress(num_chunks, num_chunks)


def unpacker(fmt, data, chunks, num_chunks):
    size = struct.calcsize(fmt)
    while size > len(data):
        chunk_num, chunk = next(chunks)
        if chunk_num % 1000 == 0:
            print_progress(chunk_num, num_chunks)
        data += chunk
    return struct.unpack(fmt, data[:size]), data[size:]


def process_read(prev_read_name, read_name, refID, pos, print_reads, all_reads):
    if all_reads or read_name[0:len(print_reads)] != print_reads:
        print(file=sys.stderr)
        hex_read_name = " ".join("{:02x}".format(ord(c)) for c in read_name)
        print(prev_read_name, read_name, refID, pos, hex_read_name, file=sys.stderr)


@contextlib.contextmanager
def smart_open(path, fallback, mode="r"):
    if path and path != "-":
        fh = open(path, mode)
    else:
        fh = fallback

    try:
        yield fh
    finally:
        if fh is not fallback:
            fh.close()


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Decompress a GZIP/BGZIP compressed file, skipping errors. Parse an uncompressed BAM.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--parse", action="store_true", help="Treat input as uncompressed BAM and parse its reads.")
    parser.add_argument("--print-reads", default="ST-", help="Print reads with names NOT matching the given string prefix to STDERR.")
    parser.add_argument("--all-reads", action="store_true", help="Print all reads to STDERR when parsing.")
    parser.add_argument("--no-cache", action="store_true", help="Do not cache decompressed chunks until the archive member is complete.")
    parser.add_argument("input_path", nargs="?", help="File to decompress or parse; '-' or none for STDIN")
    parser.add_argument("output_path", nargs="?", help="File to write decompressed results to; '-' or none for STDOUT")
    args = parser.parse_args(argv)
    return args


def chunk(input_fh):
    chunks = enumerate(iter(lambda: input_fh.read(READ_SIZE), b""), start=1)
    if input_fh is sys.stdin:
        return chunks, 0
    else:
        return chunks, int(math.ceil(os.path.getsize(input_fh.name) / READ_SIZE))


def run(args):
    with smart_open(args.input_path, sys.stdin, "rb") as input_fh, smart_open(args.output_path, sys.stdout, "wb") as output_fh:
        chunks, num_chunks = chunk(input_fh)
        if args.parse:
            parse_bam(chunks, num_chunks, args.print_reads, args.all_reads)
        else:
            for bam_chunk in optimistic_decompress(chunks, num_chunks, args.no_cache):
                output_fh.write(bam_chunk)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    run(args)
