from pycofitness.mutation.mutation import PointMutation
from .input_files_path import InputFilesPath
import unittest

from Bio import AlignIO

class PointMutationTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        #rna test files
        self.__rna_msa_file = InputFilesPath.rna_msa_file
        
        #protein test files
        self.__protein_msa_file = InputFilesPath.protein_msa_file
       
        self.__pointmutation_instance_protein = PointMutation(
            self.__protein_msa_file,
            'protein',
            num_threads=4,
            verbose=True
        )
        self.__pointmutation_instance_rna = PointMutation(
            self.__rna_msa_file,
            'rna',
            num_threads=4,
            verbose=True 
        )


    def test_delphi_epistatic_rna(self):
        """
        """
        self.__pointmutation_instance_rna.delphi_epistatic()

    def test_delphiepistatic_protein(self):
        """
        """
        self.__pointmutation_instance_protein.delphi_epistatic()


if __name__ == '__main__':
    unittest.main()
