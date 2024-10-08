from pycofitness.fasta_reader.fasta_reader import get_alignment_from_fasta_file
import ctypes 
import logging
import os
import glob 
import numpy as np


"""Python wrapper for psuedolikelihood maximization direct coupling analysis (plmDCA).
The gradient descent algorithm is implemented using c++ backend.

Authors: Mehari B. Zerihun, Fabrizio Pucci
"""

logger = logging.getLogger(__name__)


class PlmDCAException(Exception):
    """Implements exceptions related to PlmDCA computation
    """

class PlmDCA:
    """Wraps  the c++ implementation of plmDCA.

    Attributes
    ----------
        plmdca_so : str 
            Path to the plmDCA shared object created from the C++ source code
    """
    plmdca_so_paths  = glob.glob(
        os.path.join(os.path.dirname(__file__), '_plmdcaBackend*')
    )
    logger.info('\n\tplmdca backend path: {}'.format(plmdca_so_paths))
    try:
        plmdca_lib_path = plmdca_so_paths[0]
    except IndexError:
        logger.error('\n\tUnable to find plmdca dynamic library path.' 
            '\nAre you running  pydca as a module before installation?'
            '\nIn this case you need to build the plmdca shared object.'
        )
        #raise 


    def __init__(self, msa_file, biomolecule, seqid = None, lambda_h=None, 
            lambda_J = None, max_iterations = None, num_threads = None, 
            verbose = False):
        """Initilizes plmdca instances
        """
        self.__biomolecule = biomolecule.strip().upper()
        if self.__biomolecule not in ('PROTEIN', 'RNA'):
            logger.error('\n\tInvalid biomolecule type {}'.format(self.__biomolecule))
            raise PlmDCAException
        self.__msa_file = msa_file
        self.__biomolecule_int = 1 if self.__biomolecule == 'PROTEIN' else 2
        self.__num_site_states = 21 if self.__biomolecule== 'PROTEIN' else 5
        self.__num_seqs, self.__seqs_len = self._get_num_and_len_of_seqs()
        self.__seqid = 0.8 if seqid is None else seqid 
        if self.__seqid <= 0 or self.__seqid > 1.0: 
            logger.error('\n\t{} is an invalid value of sequences identity (seqid) parameter'.format(self.__seqid))
            raise PlmDCAException 
        self.__lambda_h= 10.0 if lambda_h is None else lambda_h
        if self.__lambda_h < 0 :
            logger.error('\n\tlambda_h must be a positive number. You passed lambda_h={}'.format(self.__lambda_h))
            raise PlmDCAException  
        self.__lambda_J=  6.0*(self.__seqs_len - 1) if lambda_J is None else lambda_J
        if self.__lambda_J < 0: 
            logger.error('\n\tlambda_J must be a positive number. You passed lambda_J={}'.format(self.__lambda_J))
            raise PlmDCAException
        self.__max_iterations = max_iterations if max_iterations is not None else 1000
        self.__num_threads = 1 if num_threads is None else num_threads
        self.__verbose = True if verbose else False  
        # plmdcaBackend interface
        # extern "C" float* plmdcaBackend(unsigned short const biomolecule, unsigned short const num_site_states, 
        # const char* msa_file, unsigned int const seqs_len, float const seqid, float const lambda_h, 
        # float const lambda_J, unsigned int const max_iterations)
        self.__plmdca = ctypes.CDLL(self.plmdca_lib_path)
        self.__plmdcaBackend = self.__plmdca.plmdcaBackend 
        self.__plmdcaBackend.argtypes = (ctypes.c_ushort, ctypes.c_ushort, ctypes.c_char_p, ctypes.c_uint, 
            ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_uint, ctypes.c_uint, ctypes.c_bool)
        data_size = (self.__seqs_len * (self.__seqs_len - 1) * (self.__num_site_states ** 2))/2 + self.__seqs_len * self.__num_site_states
        self.__data_size = int(data_size) 
        self.__plmdcaBackend.restype = ctypes.POINTER(ctypes.c_float * self.__data_size)
        #extern "C" void freeFieldsAndCouplings(float*& fields_and_couplings)
        self.freeFieldsAndCouplings = self.__plmdca.freeFieldsAndCouplings
        #self.freeFieldsAndCouplings.argtypes = (ctypes.POINTER(ctypes.c_float),) 
        self.freeFieldsAndCouplings.restype = None 
        log_message="""Created plmDCA instance with:
            biomolecule: {}
            MSA sequence length: {}
            total number of (unweighted) sequences: {}
            sequence identity: {}
            lambda_h: {}
            lambda_J: {}
            gradient descent iterations: {}
            number of threads: {}
        """.format(self.__biomolecule, self.__seqs_len, self.__num_seqs, 
            self.__seqid, self.__lambda_h, self.__lambda_J, self.__max_iterations,
            self.__num_threads,
        )
        logger.info(log_message)
        return None


    @property
    def biomolecule(self):
        """PlmDCA biomolecule attribute 
        """
        return self.__biomolecule


    @property
    def sequence_identity(self):
        """PlmDCA sequence_identity attribute
        """
        return self.__seqid


    @property
    def lambda_h(self):
        """PlmDCA lambda_h attribute
        """
        return self.__lambda_h
    

    @property
    def lambda_J(self):
        """PlmDCA lambda_J attribute
        """
        return self.__lambda_J


    @property
    def max_iterations(self):
        """PlmDCA max_iterations attribute
        """
        return self.__max_iterations


    @property
    def sequences_len(self):
        """PlmDCA sequences_len attribute
        """
        return self.__seqs_len


    @property
    def num_sequences(self):
        """PlmDCA num_sequences attribute
        """
        return self.__num_seqs
    
    
    @property
    def effective_num_sequences(self):
        """PlmDCA effective_num_sequences attribute
        """
        raise NotImplementedError
        

    def _get_num_and_len_of_seqs(self):
        """Obtains the length of sequences in MSA data.

        Parameters
        ----------
            self : PlmDCA 
                Instance of PlmDCA class
        
        Returns
        -------
            num_seqs, seqs_len: tuple  
                A tuple of the number and length of sequences as read from 
                the MSA file
        """
        msa_data = get_alignment_from_fasta_file(self.__msa_file)
        num_seqs = len(msa_data)
        seqs_len = len(msa_data[0])
        return num_seqs, seqs_len


    def map_index_couplings(self, i, j, a, b):
        """Couplings index mapper.

        Parameters
        ----------
            self : PlmDCA 
                An instance of PlmDCA class.
            i : int 
                Site in site-pair (i, j) such that j > i. 
            j : int 
                Site in site-pair (i, j) such that j > i.
        """
        q = self.__num_site_states
        L = self.__seqs_len
        site = int(((L *  (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * q * q)
        k =  L * q + site + b + a * q
        return k
    
    
    def get_fields_and_couplings_from_backend(self):
        """Compute fields and couplings using the C++ backend of plmDCA

        Parameters
        ----------
            self : PlmDCA
                An instance of PlmDCA class

        Returns
        -------
            fields_and_couplings : np.array
                A one-dimensional array of the fields and couplings
        """
        logger.info('\n\tComputing fields and couplings using gradient descent')
        h_J_ptr = self.__plmdcaBackend(
            self.__biomolecule_int, self.__num_site_states, self.__msa_file.encode('utf-8'), 
            self.__seqs_len,  self.__seqid, self.__lambda_h, self.__lambda_J, self.__max_iterations,
            self.__num_threads, self.__verbose
        )
        
        fields_and_couplings = np.zeros((self.__data_size,), dtype=np.float32)
        #fields_and_couplings = [c for c in h_J_ptr.contents]
        counter = 0
        for i, h_or_J in enumerate(h_J_ptr.contents):
            fields_and_couplings[i] = h_or_J
            counter += 1
        #Free fields and couplings data from PlmDCABackend
        h_J_ptr_casted = ctypes.cast(h_J_ptr, ctypes.POINTER(ctypes.c_void_p))
        self.freeFieldsAndCouplings(h_J_ptr_casted)

        #fields_and_couplings = np.array(fields_and_couplings, dtype=np.float32)
        
        try:
            assert fields_and_couplings.size == counter 
        except AssertionError:
            logger.error('\n\tData size mismatch from the plmDCA backend')
            raise
        else:
            logger.info('\n\tData size expected from plmDCA backend: {}'.format(self.__data_size))
            logger.info('\n\tData size obtained from plmDCA backend to Python: {}'.format(fields_and_couplings.size))

        return fields_and_couplings 


    def get_couplings_no_gap_state(self, fields_and_couplings_all):
        """Extract the couplings from fields and couplings numpy array

        Parameters
        ----------
            self : PlmDCA 
                An instanc of PlmDCA class
            fields_and_couplings : list/array
                A list of all fields and couplings, including gap state fields
                couplings.
        Returns
        -------
            couplings : np.array
                A one-dimensional array of the couplings
        """
        couplings = list()
        for i in range(self.__seqs_len - 1):
            for j in range(i + 1, self.__seqs_len):
                for a in range(self.__num_site_states - 1):
                    for b in range(self.__num_site_states -1):
                        indx =  self.map_index_couplings(i, j , a, b)
                        couplings.append(fields_and_couplings_all[indx])
        return np.array(couplings)

    
    def get_fields_no_gap_state(self, fields_and_couplings_all):
        """Extracts the fields from fields and couplings numpy array.

        Parameters
        ----------
            self : PlmDCA
                An instance of PlmDCA class
            fields_and_couplings : list/array
                A list of all fields and couplings, including gap state fields
                couplings.
        Returns
        --------
            fields_no_gap_state : list 
                A  list of fields excluding fields corresponding to gap states.
        """
        
        fields_all = fields_and_couplings_all[:self.__seqs_len * self.__num_site_states]
        fields_no_gap_state = list()
        for i in range(self.__seqs_len):
            for a in range(self.__num_site_states - 1): # iterate over q - 1 states to exclude gaps
                fields_no_gap_state.append(fields_all[a + i * self.__num_site_states])
        return fields_no_gap_state

    
    def get_fields_and_couplings_no_gap_state(self, fields_and_couplings_all):
        """Computes fields and couplings excluding gap state fields and gap state
        couplings.

        Parameters
        ----------
            self : PlmDCA 
                An instance of PlmDCA class
            fields_and_couplings : list/array
                A list of all fields and couplings, including gap state fields
                couplings.
        
        Returns
        --------
            fields_no_gap_state, couplings_no_gap_state : tuple
                A tuple of list of fields and couplings
        """

        
        logger.info('\n\tObtaining fields and couplings excluding gap state.')
        fields_no_gap_state = self.get_fields_no_gap_state(fields_and_couplings_all)
        couplings_no_gap_state = self.get_couplings_no_gap_state(fields_and_couplings_all)
        
        return fields_no_gap_state, couplings_no_gap_state

    
    def shift_couplings(self, couplings_ij):
        """Shifts the couplings value.

        Parameters
        ----------
            self : PlmDCA 
                An instance of PlmDCA class
            couplings_ij : np.array
                1d array of couplings for site pair (i, j)
        Returns
        -------
            shifted_couplings_ij : np.array
                A 2d array of the couplings for site pair (i, j)
        """
        qm1 = self.__num_site_states - 1
        couplings_ij = np.reshape(couplings_ij, (qm1,qm1))
        avx = np.mean(couplings_ij, axis=1)
        avx = np.reshape(avx, (qm1, 1))
        avy = np.mean(couplings_ij, axis=0)
        avy = np.reshape(avy, (1, qm1))
        av = np.mean(couplings_ij)
        couplings_ij = couplings_ij -  avx - avy + av
        return couplings_ij     


if __name__ == '__main__':
    """
    """
    
