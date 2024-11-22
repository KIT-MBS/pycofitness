from typing import Union
from Bio import AlignIO
import numpy as np 
from pycofitness.fasta_reader.fasta_reader import get_alignment_from_fasta_file
from pycofitness.resmapping.resmapping import residue_mapping
from pycofitness.fasta_reader.fasta_reader import alignment_letter2int
from logging import config
from  pycofitness.plmdca.plmdca import PlmDCA
from pycofitness.msa_trimmer.msa_trimmer import MSATrimmer
import logging 
from pathlib import Path


logger = logging.getLogger(__name__)

class SingleSiteFreqs:
    """
    """
    def __init__(self, alignment_data, seqid):
        """
        Parameters
        ----------
            alignmnet_data : np.array()
                Numpy 2d array of the alignment data, after the alignment is put in
                integer representation
            seqid : float
                Value at which beyond this sequences are considered similar. Typical
                values could be 0.7, 0.8, 0.9 and so on
        """
        self.__alignment_data = alignment_data
        seqid = 0.8 if not seqid else seqid
        self.__seqid = seqid
        
        return None 
    
    def compute_sequences_weight(self):
        """Computes weight of sequences. The weights are calculated by lumping
        together sequences whose identity is greater that a particular threshold.
        For example, if there are m similar sequences, each of them will be assigned
        a weight of 1/m. Note that the effective number of sequences is the sum of
        these weights.
        
        Returns
        -------
        
        seqs_weight : np.array()
            A 1d numpy array containing computed weights. This array has a size
            of the number of sequences in the alignment data.
        """
        logger.info('\nComputing sequence weights')
        alignment_shape = self.__alignment_data.shape
        num_seqs = alignment_shape[0]
        seqs_len = alignment_shape[1]
        seqs_weight = np.zeros((num_seqs,), dtype=np.float64)
        #count similar sequences
        for i in range(num_seqs): #parallel_range(num_seqs):
            seq_i = self.__alignment_data[i]
            for j in range(num_seqs):
                seq_j = self.__alignment_data[j]
                iid = np.sum(seq_i==seq_j)
                if np.float64(iid)/np.float64(seqs_len) > self.__seqid:
                    seqs_weight[i] += 1
        #compute the weight of each sequence in the alignment
        for i in range(num_seqs): seqs_weight[i] = 1.0/float(seqs_weight[i])
        return seqs_weight

    def compute_single_site_freqs(self, num_site_states):
        """Computes single site frequency counts for a particular aligmnet data.

        Parameters
        ----------
            num_site_states : int
                An integer value fo the number of states a sequence site can have
                including a gap state. Typical value is 5 for RNAs and 21 for
                proteins.

        Returns
        -------
            single_site_freqs : np.array()
                A 2d numpy array of of data type float64. The shape of this array is
                (seqs_len, num_site_states) where seqs_len is the length of sequences
                in the alignment data.
        """
        alignment_shape = self.__alignment_data.shape
        seqs_weight = self.compute_sequences_weight()
        logger.info('\nComputing single site frequencies (normalized by Meff)')
        #num_seqs = alignment_shape[0]
        seqs_len = alignment_shape[1]
        m_eff = np.sum(seqs_weight)
        single_site_freqs = np.zeros(shape = (seqs_len, num_site_states),
            dtype = np.float64)
        for i in range(seqs_len):
            for a in range(1, num_site_states + 1):#we need gap states single site freqs too
                column_i = self.__alignment_data[:,i]
                freq_ia = np.sum((column_i==a)*seqs_weight)
                single_site_freqs[i, a-1] = freq_ia/m_eff
        return single_site_freqs
    
    def get_reg_single_site_freqs(self, seqs_len,num_site_states, pseudocount = None):
        """Regularizes single site frequencies.

        Parameters
        ----------
            seqs_len : int
                The length of sequences in the alignment data
            num_site_states : int
                Total number of states that a site in a sequence can accommodate. It
                includes gap states.
            pseudocount : float
                This is the value of the relative pseudo count of type float.
                theta = lambda/(meff + lambda), where meff is the effective number of
                sequences and lambda is the real pseudo count.

        Returns
        -------
            reg_single_site_freqs : np.array()
                A 2d numpy array of shape (seqs_len, num_site_states) of single site
                frequencies after they are regularized.
        """
        logger.info('\nRegularizing single site frequencies')
        if pseudocount is None: pseudocount = 0.5
        if pseudocount < 0.0 or pseudocount > 1.0:
            logger.error(f'\nIvalid pseudocount value {pseudocount}. Must be between zero and one')
            raise ValueError
        single_site_freqs = self.compute_single_site_freqs(num_site_states=num_site_states)
        reg_single_site_freqs = single_site_freqs
        theta_by_q = np.float64(pseudocount)/np.float64(num_site_states)
        for i in range(seqs_len):
            for a in range(num_site_states):
                reg_single_site_freqs[i, a] = theta_by_q + (1.0 - pseudocount)*reg_single_site_freqs[i, a]
        return reg_single_site_freqs


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
        
        self.__refseq = None # will be set to refseq at after methods call
        
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
    
    
    def delphi_epistatic(self, save_path: Union[None, str]=None):
        """change in energy function is computed from fields and couplings
        """ 

        # Init container for h and J coefficients
        hJ_coeff_lines = []

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
                    hJ_coeff_lines.append(f"{b} {i+1} {j+1} {delta_hi_ba} {delta_Jij}") # add +1 to residue-ids to match the fasta format conventions
                for j in range(i + 1, self.__num_sites): # j > i sites
                    res_at_j = self.__refseq[j]
                    if res_at_j not in list_of_standard_residues: continue 
                    couplings_loc_ij = int(((L *  (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * self.__qm1 * self.__qm1)
                    k_mutant = couplings_loc_ij + res_mapping[res_at_j]  + res_mapping[b] * self.__qm1
                    k_background = couplings_loc_ij + res_mapping[res_at_j] + res_mapping[res_background] * self.__qm1
                    delta_Jij = self.__couplings[k_mutant] - self.__couplings[k_background]
                    delta_Jij_total += delta_Jij
                    hJ_coeff_lines.append(f"{b} {i+1} {j+1} {delta_hi_ba} {delta_Jij}") # add +1 to residue-ids to match the fasta format conventions
                deltaPhi = delta_Jij_total + delta_hi_ba
                mut_data[i][b] = deltaPhi

        # Save file of h and J coefficients if required
        if save_path is not None:
            with open(save_path, "w") as fs:
                fs.write("\n".join(hJ_coeff_lines))

        return mut_data                         


    def delphi_independent_site(self, pseudocount=None, save_path: Union[None, str]=None):
        """ 
        """

        # Init container for h coefficients
        h_coeff_lines = []

        alignment_data = alignment_letter2int(self.__msa_file, biomolecule=self.__biomolecule)
        res2int = residue_mapping[self.__biomolecule]
        res2lett = {val: key for (key, val) in res2int.items()}
        refseq_int = alignment_data[0]
        seqs_len = len(refseq_int)
        refseq_str = ''.join([res2lett[res] for res in refseq_int])
        self.__refseq = refseq_str
        self.__num_sites  = len(self.__refseq)
        num_site_states = 5 if self.__biomolecule == 'RNA' else 21
        reg_single_site_freqs = SingleSiteFreqs(
            alignment_data=alignment_data, 
            seqid=self.__seqid
        ).get_reg_single_site_freqs(num_site_states= num_site_states, pseudocount=pseudocount, seqs_len=seqs_len)
        
        logger.info(f'\nPerforming point mutation using independent site model for sequence {refseq_str}')
        fields_vec = np.log(reg_single_site_freqs)
        
        deltas_dict = dict()
        for i in range(seqs_len):
            deltas_dict[i] = dict()
            a = self.refseq[i]
            ref_allele_indx = res2int[a]
            for mut_allele_indx in range(num_site_states -1):
                b = res2lett[mut_allele_indx]
                curr_delta = fields_vec[i, mut_allele_indx] - fields_vec[i, ref_allele_indx]
                deltas_dict[i][b] = curr_delta
                h_coeff_lines.append(f"{b} {i+1} {curr_delta}") # add +1 to residue-ids to match the fasta format conventions

        # Save file of h coefficients if required
        if save_path is not None:
            with open(save_path, "w") as fs:
                fs.write("\n".join(h_coeff_lines))

        return deltas_dict
                
                

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
    test_msa_RNA = '/home/mehari/pycofitness_github_cloned/pycofitness/tests/tests_input/trimmed_50RNAseqs.fa'
    biomolecule= 'rna' 
       
    point_mut_inst = PointMutation(test_msa_RNA, biomolecule, max_iterations=20000, num_threads=6, verbose=True)
    #delphi_epi=point_mut_inst.delphi_epistatic()
    delphi_ind = point_mut_inst.delphi_independent_site()
    delphi = delphi_ind
    counter = 0
    for i in range(len(point_mut_inst.refseq)):
        for a in point_mut_inst.standard_residues:
            counter += 1
            print(i + 1, delphi[i][a])
            
    print(point_mut_inst.standard_residues)
    
    return None

def main():
    """ 
    """
    logger.info(f'\n\tYou are running {__file__}')
    test_rna()
    
    
                        
if __name__ == '__main__':
    from pycofitness.main import configure_logging
    configure_logging()
    test_rna()
