#!/usr/bin/env -S uv run --quiet

# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "biopython",
# ]
# ///
from Bio.Seq import Seq

my_seq = Seq("GATTACA")
print(my_seq)

rc = my_seq.reverse_complement()
print(rc)
