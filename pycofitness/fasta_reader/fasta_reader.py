from Bio import AlignIO
import logging 

logger = logging.getLogger(__name__)

def get_alignment_from_fasta_file(file_name):
    """Read sequences from FASTA file using Bio.AlignIO.read()
    Parameters
    ----------
        file_name : str
            Path to FASTA formatted file.
    Returns
    -------
        alignment : list
            A list of biomolecular sequence strings.
    """
    alignment = []
    try:
        record_iterator = AlignIO.read(file_name, 'fasta')
        #biopython just reads the records if there are tags (>some key).
        #It doesn't know if the file is really a biological sequence or not
    except Exception as expt:
        error_msg='\n\tError occured while reading from fasta file: {}.' +\
            '\n\tError type:{}\n\tArguments:{!r}'
        logger.error(error_msg.format(file_name, type(expt).__name__, expt.args))
        raise
    else:
        if any(True for _ in record_iterator):
            for record in record_iterator:
                seq = record.seq.strip()
                if seq: alignment.append(seq.upper())
            if not alignment:
                logger.error(
                    '\n\trecord_iterator returned by AlignIO.read()'
                    ' has no sequences',
                )
                raise ValueError

        else:
            logger.error(
                '\n\trecord_iterator returned by AlignIO.read() is empty',
            )
            raise ValueError
    return alignment