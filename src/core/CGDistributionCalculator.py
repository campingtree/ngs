from collections import Counter
from math import floor
from typing import Dict, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class CGDistributionCalculator:
    def __init__(self, data_path: str):
        self.__seq = SeqIO.parse(data_path, 'fastq')

    def calculate(self) -> Dict[int, List[SeqRecord]]:
        frequencies = dict.fromkeys(range(101), [])

        for record in self.__seq:
            seq_len = len(record.seq)

            if seq_len == 0:
                continue

            chars = dict(Counter(record))
            if all(c in chars for c in ('C', 'G')):
                freq = floor( (chars['C'] + chars['G']) * 100 / seq_len)
                frequencies[freq] = frequencies[freq] + [record]

        return frequencies
