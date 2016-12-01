#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import sys
import traceback
import zlib
import os
import math
import contextlib


READ_SIZE = 4096
DECOMPRESSOR_OPTIONS = zlib.MAX_WBITS | 16  # max window, gzip header/trailer
GZIP_MAGIC = "\037\213"


def optimistic_decompress(input_path, output_path):
    with open(input_path, "rb") as in_fh, smart_open(output_path, sys.stdout, "wb") as out_fh:
        num_chunks = int(math.ceil(os.path.getsize(input_path) / READ_SIZE))
        decompressor = zlib.decompressobj(DECOMPRESSOR_OPTIONS)
        unused = b""
        chunks = enumerate(iter(lambda: in_fh.read(READ_SIZE), b""), start=1)
        for chunk_num, chunk in chunks:
            try:
                bam_chunk = decompressor.decompress(unused + chunk)
                out_fh.write(bam_chunk)
            except Exception:
                print_exception(chunk_num, decompressor)
                unused = find_next_member(chunks, decompressor)
                decompressor = zlib.decompressobj(DECOMPRESSOR_OPTIONS)
                continue

            if decompressor.unconsumed_tail:
                print("unexpected unconsumed tail of length {}".format(len(decompressor.unconsumed_tail)), file=sys.stderr)

            unused = decompressor.unused_data
            if unused:
                decompressor = zlib.decompressobj(DECOMPRESSOR_OPTIONS)

            if chunk_num % 1000 == 0:
                print_progress(chunk_num, num_chunks)

        # gzip module only does this at the end, not on new members
        bam_chunk = decompressor.flush()
        out_fh.write(bam_chunk)

        print_progress(num_chunks, num_chunks)


def print_exception(chunk_num, decompressor):
    error_offset = chunk_num * READ_SIZE - len(decompressor.unconsumed_tail)
    print(file=sys.stderr)
    print("decompression error at offset {} (unused data: {})".format(error_offset, len(decompressor.unused_data)), file=sys.stderr)
    traceback.print_exc(file=sys.stderr)


def find_next_member(chunks, decompressor):
    next_part = decompressor.unconsumed_tail.find(GZIP_MAGIC)
    unused = b""
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
    percent_complete = chunk_num / num_chunks
    print("Decompressed chunk {:{width}}/{:{width}} ({:.1%})".format(
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


if __name__ == "__main__":
    optimistic_decompress(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
