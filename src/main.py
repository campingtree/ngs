import matplotlib.pyplot as plt

from core import CGDistributionCalculator, FASTQQualityEncodingFormatDetector
from utils import BlastUtils

if __name__ == '__main__':
    data_path = '../data/reads.fastq'

    # a
    quality_enc_format = FASTQQualityEncodingFormatDetector(data_path).detect()
    print(f'Detected quality encoding format:  {quality_enc_format}')

    # b
    frequencies = CGDistributionCalculator(data_path).calculate()
    print('Outputting C/G distribution graph to a file...')
    x_values = list(frequencies.keys())
    y_values = [len(seq) for seq in frequencies.values()]
    plt.plot(x_values, y_values)
    plt.title('C/G distribution in reads')
    plt.xlabel('% of C/G nucleotides')
    plt.ylabel('Number of reads')
    plt.grid()
    plt.savefig('cg_distribution.png')

    # c
    print('Doing BLAST search for 5 sequences from each specific peak...')
    peak_records = []
    peak_record_selection_criteria = [
        (31, 1), (33, 2), (35, 2),          # peak 1
        (50, 1), (52, 2), (54, 1), (56, 1), # peak 2
        (66, 1), (68, 1), (70, 2), (72, 1)  # peak 3
    ]

    for freq, num in peak_record_selection_criteria:
        peak_records += frequencies[freq][:num]

    with open('blast_results.txt', 'w') as f:
        for peak_rec in peak_records:
            micro_org_name = BlastUtils.blast_search_first_match(peak_rec.seq)
            f.write(f'{peak_rec.name}\n{micro_org_name}\n\n')
