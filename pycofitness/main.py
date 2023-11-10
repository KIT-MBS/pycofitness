from pycofitness.logging_config.logging_config import LOGGING_CONFIG
from pycofitness.logging_config.logging_config import ConsoleColor as c_color
from pycofitness.mutation import mutation
from pycofitness.mutation import mutation 
from argparse import ArgumentParser
import pathlib
import logging 
import sys 


logger = logging.getLogger(__name__)

class CmdArgs:
    """Defines command line argument variables for plmDCA and pycofitness.
    """
    msa_file = 'msa_file'
    msa_file_help = """Multiple sequence alignment (MSA) file in FASTA format.
    """
    biomolecule = 'biomolecule'
    biomolecule_help = """Type of biomolecule. It should be either protein or RNA 
    in lower or upper case letters.
    """
    verbose_optional = '--verbose'
    verbose_optional_help = """Show logging information on the terminal.
    """
    seqid_optional = '--seqid'
    seqid_optional_help = """Cut-off value of sequences similarity above which they
    are lumped together.
    """
    lambda_h_optional = '--lambda_h'
    lambda_h_optional_help = """Value of fields penalizing constant for L2 
    regularization of fields.
    """
    lambda_J_optional = '--lambda_J'
    lambda_J_optional_help = """Value of couplings penalizing constant for L2 
    regularization of couplings.
    """
    max_iterations_optional = '--max_iterations'
    max_iterations_help = """Number of iterations for gradient descent 
    for negative pseudolikelihood minimization.
    """
    num_threads_optional = '--num_threads'
    num_threads_help = "Number of threads from plmDCA computation"
    output_dir_optional = '--output_dir'
    output_dir_help = """Directory path to which output results are written.
    If the directory is not existing, it will be created. If this path is not provided, an output directory
    is created using the base name of the MSA file, with "output_" prefix
    added to it.
    """     
#End of class CmdArgs

logger = logging.getLogger(__name__)

def configure_logging():
    """Configures logging. When configured, the logging level is INFO and
    messages are logged to stream handler. Log level name are colored whenever
    the terminal supports that. INFO level is Green, WARNING level is Yellow and
    ERROR level is Red.
    """
    logging.config.dictConfig(LOGGING_CONFIG)
    logging.addLevelName(logging.INFO, '{}{}{}'.format(
        c_color.green, logging.getLevelName(logging.INFO), c_color.nocolor))
    logging.addLevelName(logging.WARNING, '{}{}{}'.format(
        c_color.yellow, logging.getLevelName(logging.WARNING), c_color.nocolor))
    logging.addLevelName(logging.ERROR, '{}{}{}'.format(
        c_color.red, logging.getLevelName(logging.ERROR), c_color.nocolor))
    return None

    
def execute_from_command_line(biomolecule, msa_file, 
        seqid = None, lambda_h = None, lambda_J = None, 
        max_iterations = None ,verbose = False, output_dir = None,
        num_threads = None):
    
    if verbose : configure_logging()
    
    mut_inst = mutation.PointMutation(msa_file, biomolecule, seqid=seqid, lambda_h=lambda_h, 
        lambda_J=lambda_J, num_threads=num_threads,
        max_iterations=max_iterations, verbose=verbose
    )
    delphi_dict= mut_inst.delphi_epistatic()
    if output_dir:
        dir_output = pathlib.Path(output_dir)
        dir_output.mkdir()
    else:
        msa_file_base = pathlib.Path(msa_file).stem 
        dir_output = pathlib.Path(f'output_{msa_file_base}')
        dir_output.mkdir()
    outfile = dir_output / ('deltaphi.txt')
    with open(outfile ,'w') as fh:
        fh.write('#site\treference\talternative\tscore\n')
        for i in range(mut_inst.num_sites):
            for res in mut_inst.standard_residues:
                fh.write('{}\t{}\t{}\t{}\n'.format(i+1,mut_inst.refseq.upper()[i], res, delphi_dict[i][res]))
    return None 


def run_mutation():
    """
    """
    parser = ArgumentParser()
    
    parser.add_argument(CmdArgs.biomolecule, help=CmdArgs.biomolecule_help)
    parser.add_argument(CmdArgs.msa_file, help=CmdArgs.msa_file_help)
    parser.add_argument(CmdArgs.seqid_optional, help=CmdArgs.seqid_optional_help, type=float)
    parser.add_argument(CmdArgs.lambda_h_optional, help=CmdArgs.lambda_h_optional_help, type=float)
    parser.add_argument(CmdArgs.lambda_J_optional, help=CmdArgs.lambda_J_optional_help, type=float)
    parser.add_argument(CmdArgs.max_iterations_optional, help=CmdArgs.max_iterations_help, type=int)
    parser.add_argument(CmdArgs.num_threads_optional, help=CmdArgs.num_threads_help, type=int)
    parser.add_argument(CmdArgs.verbose_optional, help=CmdArgs.verbose_optional_help, action='store_true')
    parser.add_argument(CmdArgs.output_dir_optional, help=CmdArgs.output_dir_help)   
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    args_dict = vars(args)
    
    execute_from_command_line(
        args_dict.get('biomolecule'), 
        args_dict.get('msa_file'), 
        seqid=args_dict.get('seqid'),
        lambda_h=args_dict.get('lambda_h'),
        lambda_J = args_dict.get('lambda_J'),
        max_iterations = args_dict.get('max_iterations'),
        num_threads = args_dict.get('num_threads'),
        output_dir = args_dict.get('output_dir'),
        verbose = args_dict.get('verbose'),
    )

    
    return None 


if __name__ == '__main__':
    run_mutation()


