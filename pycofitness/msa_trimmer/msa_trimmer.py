from Bio import AlignIO
from pycofitness.resmapping.resmapping import residue_mapping
import logging

"""Removes gaps from MSA data based of the reference sequence. 
"""

logger = logging.getLogger(__name__)


class MSATrimmer:
    """
    """
    def __init__(self, msa_file, biomolecule):
        """
        """
        self.__msa_file = msa_file 
        self.__biomolecule = biomolecule
        self.__alignment_data = self.get_alignment_data()
        return None 
                
    
    def get_alignment_data(self):
        alignment_data = list(AlignIO.read(self.__msa_file, 'fasta'))
        seq_data_list = list()
        for record in alignment_data:
            seqid = record.id
            seq  = record.seq.upper()
            seq_data_list.append((seqid, str(seq)))
        return seq_data_list
    
    def get_trimmed_cols(self):
        """
        """
        mapping = residue_mapping[self.__biomolecule]
        std_residues = list(mapping.keys())
        trimmed_cols = list()
        refseq = self.__alignment_data[0][1]
        for i, res in enumerate(refseq):
            if res in std_residues: trimmed_cols.append(i)
        return trimmed_cols
       
    def get_trimmed_msa_data(self):
        """Removes all non-standard residues
        """
        trimmed_cols = self.get_trimmed_cols()
        trimmed_msa_data = list()
        for seq_id, seq in self.__alignment_data:
            trimmed_seq = list()
            for i in trimmed_cols:
                trimmed_seq.append(seq[i])
            trimmed_seq = ''.join(trimmed_seq)
            trimmed_msa_data.append((seq_id, trimmed_seq))
        return trimmed_msa_data
    

def main():
    """ 
    """
    test_seq_RNA = '../../tests/tests_input/test_trim_RNA.fa'
    msa_trimmer = MSATrimmer(msa_file=test_seq_RNA, biomolecule='RNA')
    msa_trimmer.get_trimmed_msa_data()
    return None 

if __name__ == '__main__':
    """
    """
    main()
