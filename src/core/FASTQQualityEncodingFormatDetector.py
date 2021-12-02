class FASTQQualityEncodingFormatDetector:
    def __init__(self, data_path: str):
        self.__ascii_records = []
        with open(data_path, 'rU') as f:
            for line in f:
                next(f)
                next(f)
                qual_ascii = next(f).rstrip()
                ascii_list = list(map(ord, qual_ascii))
                self.__ascii_records.append(ascii_list)

    def detect(self) -> str:
        checks = []

        for record in self.__ascii_records:
            min_q = min(record)
            max_q = max(record)
            checks.append(min_q)
            checks.append(max_q)

        min_q = min(checks)
        max_q = max(checks)

        if 64 > min_q >= 33 and max_q == 74:
            return 'Illumina 1.8+ Phred+33'
        elif min_q >= 64 and 74 < max_q <= 104:
            return 'Illumina 1.3+ Phred+64'
        elif 64 > min_q >= 33 and max_q <= 73:
            return 'Sanger Phred+33'
        elif 73 > min_q >= 59 and max_q <= 104:
            return 'Solexa Solexa+64'
        elif min_q >= 67 and 73 < max_q <= 105:
            return 'Illumina 1.5+ Phred+64'
        else:
            return 'Quality encoding format not recognised'
