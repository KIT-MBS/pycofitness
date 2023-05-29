from Bio import AlignIO
from pycofitness.fasta_reader.fasta_reader import get_alignment_from_fasta_file
from pycofitness.resmapping.resmapping import residue_mapping
from logging import config
from  pycofitness.plmdca.plmdca import PlmDCA
from pycofitness.msa_trimmer.msa_trimmer import MSATrimmer
import logging 
from pathlib import Path


logger = logging.getLogger(__name__)


class PlmDCAParams:
    """
    """
    def __init__(self, msa_file:str, biomolecule:str) -> None:
        self.__msa_file = msa_file
        self.__biomolecule = biomolecule
        self.__refseq = get_alignment_from_fasta_file(self.__msa_file)[0]
        return None 
    
    @property
    def refseq(self):
        """
        """
        return self.__refseq
    
    
    def plmdca_params(self, seqid:float=None, lambda_h:float=None, lambda_J:float=None,
            max_iterations:int=None, num_threads:int=None, verbose:bool =False)->tuple:
        """
        """
        logger.info('\n\tComputing fields and couplings using plmDCA algorithm')
        plmDCA_inst =  PlmDCA(self.__msa_file, self.__biomolecule, seqid=seqid, 
            lambda_h=lambda_h, lambda_J=lambda_J, max_iterations=max_iterations,
            num_threads=num_threads, verbose=verbose,
        )
        
        
        fields_and_couplings_all = plmDCA_inst.get_fields_and_couplings_from_backend()
        return plmDCA_inst.get_fields_and_couplings_no_gap_state(fields_and_couplings_all)
    

class PointMutation:
    """
    """
    def __init__(self, msa_file, biomolecule, seqid=None, lambda_h=None, lambda_J=None, num_threads=None,
                 max_iterations=None, verbose=False):
        """ 
        """
        self.__biomolecule = biomolecule.strip().upper() 
        self.__msa_file = self.get_trimmed_msa_file(msa_file, self.__biomolecule)
        self.__seqid = seqid
        self.__lambda_h = lambda_h
        self.__lambda_J = lambda_J
        self.__num_threads = num_threads
        self.__max_iterations = max_iterations
        self.__verbose = verbose
        self.__qm1 = 4 if self.__biomolecule == 'RNA' else 20
        
        return None 
    @property
    def num_sites(self):
        return self.__num_sites
    
    @property
    def refseq(self):
        return self.__refseq 
    
    @property
    def standard_residues(self):
        """ 
        """
        return list(residue_mapping[self.__biomolecule].keys())
    
    def get_trimmed_msa_file(self, orig_msa_file, biomolecule):
        """ 
        """
        alignmed_data = MSATrimmer(msa_file=orig_msa_file, biomolecule=biomolecule).get_trimmed_msa_data()
        trimed_msa_file =  'trimmed_' + Path(orig_msa_file).name
        with open(trimed_msa_file, 'w') as fh:
            for sid, seq in alignmed_data:
                fh.write('>{}\n{}\n'.format(sid, seq))
        return trimed_msa_file
    
    def initialize_fields_and_couplings(self):
        """
        """
        __params_inst= PlmDCAParams(self.__msa_file, self.__biomolecule)
        self.__fields, self.__couplings = __params_inst.plmdca_params(
            seqid=self.__seqid, lambda_h=self.__lambda_h, lambda_J=self.__lambda_J,
            max_iterations=self.__max_iterations,num_threads=self.__num_threads, verbose=self.__verbose
        )
        self.__refseq = __params_inst.refseq
        self.__num_sites = len(self.__refseq)
        return None 
    

    def delphi_epistatic(self):
        """change in energy function is computed from fields and couplings
        """ 
        self.initialize_fields_and_couplings()
        logger.info('\n\tPerforming epistatic model point mutation for sequence\n\t{}'.format(self.__refseq))
        res_mapping =  residue_mapping[self.__biomolecule]
        list_of_standard_residues = list(res_mapping.keys())
        mut_data = dict()
        L = self.__num_sites
        for i in range(self.__num_sites):
            res_background = self.__refseq[i] 
            if res_background not in list_of_standard_residues:
                logger.error(f'\nReference sequence contains non-standard residue at position {i}')
                raise ValueError
            mut_data[i] = dict()
            fields_loc_i = i * self.__qm1
            hi_background = self.__fields[fields_loc_i + res_mapping[res_background]]
            for b in list_of_standard_residues: # carry out mutations at site i
                hi_mutant = self.__fields[fields_loc_i + res_mapping[b]]
                delta_hi_ba = hi_mutant - hi_background
                delta_Jij_total = 0
                for j in range(i): # j < i sites
                    res_at_j = self.__refseq[j]
                    if res_at_j not in list_of_standard_residues: continue
                    couplings_loc_ij = int(((L *  (L - 1)/2) - (L - j) * ((L-j)-1)/2  + i  - j - 1) * self.__qm1 * self.__qm1)
                    k_mutant = couplings_loc_ij + res_mapping[b] + res_mapping[res_at_j] * self.__qm1
                    k_background = couplings_loc_ij + res_mapping[res_background] + res_mapping[res_at_j] * self.__qm1 
                    delta_Jij = self.__couplings[k_mutant] - self.__couplings[k_background]
                    delta_Jij_total += delta_Jij
                for j in range(i + 1, self.__num_sites): # j > i sites
                    res_at_j = self.__refseq[j]
                    if res_at_j not in list_of_standard_residues: continue 
                    couplings_loc_ij = int(((L *  (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * self.__qm1 * self.__qm1)
                    k_mutant = couplings_loc_ij + res_mapping[res_at_j]  + res_mapping[b] * self.__qm1
                    k_background = couplings_loc_ij + res_mapping[res_at_j] + res_mapping[res_background] * self.__qm1
                    delta_Jij = self.__couplings[k_mutant] - self.__couplings[k_background]
                    delta_Jij_total += delta_Jij
                deltaPhi = delta_Jij_total + delta_hi_ba
                mut_data[i][b] = deltaPhi
        return mut_data                         

#End of class PointMutation
def test_protein():
    """
    """
    test_msa_Protein = 'TestMSA/TPK1.faclean'
    biomolecule= 'protein'  # or 'rna'
       
    point_mut_inst = PointMutation(test_msa_Protein, biomolecule, max_iterations=20000, num_threads=6, verbose=True)
    delphi_epi=point_mut_inst.delphi_epistatic()
    #delphi_ind = point_mut_inst.delphi_independent()
    counter = 0
    for i in range(len(point_mut_inst.refseq)):
        for a in point_mut_inst.standard_residues:
            counter += 1
            print(counter, delphi_epi[i][a])
    return None

def test_rna():
    test_msa_Protein = 'TestMSA/RF00059.faclean'
    biomolecule= 'rna' 
       
    point_mut_inst = PointMutation(test_msa_Protein, biomolecule, max_iterations=20000, num_threads=6, verbose=True)
    delphi_epi=point_mut_inst.delphi_epistatic()
    #delphi_ind = point_mut_inst.delphi_independent()
    counter = 0
    for i in range(len(point_mut_inst.refseq)):
        for a in point_mut_inst.standard_residues:
            counter += 1
            print(counter, delphi_epi[i][a])
            
    print(point_mut_inst.standard_residues)
    
    return None

def main():
    """ 
    """
    logger.info(f'\n\tYou are running {__file__}')
    test_rna()
    
    
                        
if __name__ == '__main__':
    from pycofitness.logging_config.logging_config import configure_logging
    configure_logging()
    main()
