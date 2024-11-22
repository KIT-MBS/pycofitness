<img width="322" alt="222" src="Logo.png">



# About `pycofitness`
`pycofitness` is a Python implementation for *in silico* mutagenesis of a target protein or RNA sequence, given as input a multiple sequence alignment (MSA) of    homologs to the target. It employs a **__pseudolikelihood maximization__**  direct coupling analysis (DCA) algorithm to infer a coevolutionary energy and then uses it to predict the impact of all single-site mutations in the target. 

The input file for `pycofitness` must be an MSA  in FASTA format, with the first sequence representing the target sequence to mutate. The output consists of a txt file reporting the change of fitness for all single-site variants of the target. 

More information about pycofitness can be found in the [associated paper](https://academic.oup.com/bioinformatics/article/40/2/btae074/7604531). When using pycofitness, please cite F. Pucci, M. Zerihun, M. Rooman, A. Schug, *"pycofitness---Evaluating the fitness landscape of RNA and protein sequences"*, Bioinformatics 40, btae074 (2024).

# Prerequisites
`pycofitness` is implemented mainly in Python, with the pseudolikelihood maximization parameter inference part implemented using C++ backend to enable parallelization. It requires: 
* Python 3
* C++ compiler that supports C++11, e.g., the GNU compiler collection (GCC)
* Optionally, OpenMP for multithreading support

# Installing

`pycofitness` can be installed from PyPI or from this repository. To install the current version of `pycofitness` from PyPI, run: 

```bash
$ pip install pycofitness
```
The alternative way is to clone this repository and install pycofitness as:

```bash
$ git clone https://github.com/KIT-MBS/pycofitness/
$ cd pycofitness
$ python setup.py install
```


# Running `pycofitness` from the command line
When `pycofitness` is installed, it provides the `pycofitness` command  that can be executed from the command line.
For example:
```bash
$ pycofitness <biomolecule> <msa_file> [--optional_arguments]
``` 
where `<biomolecule>` can take the values "protein" or "RNA" (case insensitive), and `<msa_file>` is the MSA file in  fasta format. 
The optional arguments of pycofitness input are:

`-h, --help:` Shows  command line information.

`--seqid SEQID:` Cut-off value of sequence similarity, above which sequences are lumped together in the calculation of frequencies.

`--lambda_h LAMBDA_h:` Value of the penalizing constant for L2 regularization of fields in the pseudolikelyhood maximization DCA inference step (see [here](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/36/7/10.1093_bioinformatics_btz892/2/bioinformatics_36_7_2264_s2.pdf?Expires=1702572670&Signature=02fkMyK1WmMFw69v-CfRjpNnzeLsLetV7xNIyi6RGIbgMYTyWjckjd4jxtF6XseVwe5E8JL2v4mWdUXm26C5pMtl5zlaN8zrWDanolXkgLI6dMK~9DvP-mZtEbQus49g34~wi7w~nXbBBtdzzlyFYLTlM1HIMn8i2TRzVAEKECdq~4UAccxZ1MIo1-A-fhsBqb8ZS0n7wqeimPFimgq~Tvi3nmiI1h0ud7eNh7JSaDQ-WPdIKRACOPEd1m1w5EP79NqgUuSlQvuKxnHvORaWwdTcZW0EtLpYk5-TtJWxU5szujvlrFCnSeFDeDWpX5darWr~O8Q35NfZaUsi0N8yCw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA) for  details).

`--lambda_J LAMBDA_J:` Value of the penalizing constant for L2 regularization of couplings in the pseudolikelyhood maximization DCA inference step (see [here](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/36/7/10.1093_bioinformatics_btz892/2/bioinformatics_36_7_2264_s2.pdf?Expires=1702572670&Signature=02fkMyK1WmMFw69v-CfRjpNnzeLsLetV7xNIyi6RGIbgMYTyWjckjd4jxtF6XseVwe5E8JL2v4mWdUXm26C5pMtl5zlaN8zrWDanolXkgLI6dMK~9DvP-mZtEbQus49g34~wi7w~nXbBBtdzzlyFYLTlM1HIMn8i2TRzVAEKECdq~4UAccxZ1MIo1-A-fhsBqb8ZS0n7wqeimPFimgq~Tvi3nmiI1h0ud7eNh7JSaDQ-WPdIKRACOPEd1m1w5EP79NqgUuSlQvuKxnHvORaWwdTcZW0EtLpYk5-TtJWxU5szujvlrFCnSeFDeDWpX5darWr~O8Q35NfZaUsi0N8yCw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA) for  details).

`--max_iterations MAX_ITERATIONS:` Maximum number of iterations for the gradient descent in the negative pseudolikelihood minimization step.

 `--num_threads NUM_THREADS:` Number of threads used in the computation.

 `--output_dir OUTPUT_DIR:` Directory path to which output results are written. If the directory is not existing, it is created. If this path is not provided, an output directory is created using the base name of the MSA file, with the "output_" prefix added to it.

  `--hJ_path hJ_PATH:` File path to save h and J coefficients (optional, they will not be saved by default).

 `--verbose:` Show logging information on the terminal.
 

# Using `pycofitness` as a Python library
After installation, `pycofitness` can be imported into other Python source codes, and used. For example,  

```python 
from pycofitness.mutation import PointMutation

point_mutation = PointMutation(msa_file, biomolecule)
deltas = point_mutation.delphi_epistatic()
#Print the values of delta_dict by iteration over sites and bases/residues

for i in range(len(deltas)):
    for res in deltas[i]:
        print(i + 1, res, deltas[i][res])
#Note: we added 1 to i to count sites starting from one instead of zero.
```
It is also possble to pass parameters (e.g., the number of iterations or the number of threads for parallel execution)  as 
keyword arguments to the `PointMutation` constructor:
```python 
point_mutation = PointMutation(msa_file, biomolecule,
    max_iterations = 1000,
    num_threads = 4,
    seqid = 0.9,
    lambda_J = 10.0,
    lambda_h = 5.0,
    hJ_path = './hJ_coeff.txt',
    verbose = True
)
```
where `max_iterations` is the number of maximum iterations for gradient descent, `num_threads` is the number of 
threads for parallel execution, `seqid` is the sequence identity threshold value (if sequences have a higher similarity  than this value, they are regarded as identical), `lambda_J` and `lambda_h` are penalizing constants for L2 regularization, `hJ_path` is the file path where the coefficients h and J will be saved. If `verbose` is set to `True,` logging information is given.

# Preprocessing input MSAs

The input data for pycofitness is an MSA file in FASTA format. We chose not to include tools to generate MSA in pycofitness because, in this way, users can employ
their favorite MSA alignment and curation methods. There is only one compulsory adjustment within the input MSA: the sequence
to mutate (the target sequence) has to be positioned as the first one in the MSA. Note that pycofitness automatically removes MSA columns containing non-standard residues and gaps with respect to the target sequence.

# Results interpretation

Here is a few-line example of a pycofitness output file:

#site	reference	alternative	score<br> 1	M	A	-0.14<br>1	M	C	-3.54<br>1	M	D	-5.12<br>1	M	E	-7.64<br>....<br>....<br>

The first column is the position in the target sequence, the second column is the wild-type amino acid or nucleic acid base, and the third column is the substituted amino acid or nucleic acid base. In the last column, the pycofitness score is provided. Negative values of this score correspond to variants that are less fit than the wild type, and  positive values mean variants that  are  fitter than wild type.

# Version update

Modifications in version 1.4: The default regularization parameters have been changed to lambda_J = 6.0*(L-1) and lambda_h = 10.0.  


