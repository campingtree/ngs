from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord


class BlastUtils:

    @staticmethod
    def blast_search_first_match(seq: SeqRecord) -> str:
        result_handle = NCBIWWW.qblast('blastn', 'nt', seq,
            alignments=1,
            hitlist_size=1,
            entrez_query='txid2[ORGN]') # bacteria
        records = NCBIXML.parse(result_handle)
        for rec in records:
            for align in rec.alignments:
                return align.hit_def
