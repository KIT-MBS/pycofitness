<img width="322" alt="222" src="Logo.png">



# About `pycofitness`
`pycofitness` is a Python implementation of *in silico* mutagenesis for multiple sequence alignments (MSA) of protein or RNA sequence families using **__pseudolikelihood maximization__**  direct couplings analysis (DCA) algorithm. The input file for `pycofitness` must be an MSA of sequences in fasta format, with the first sequence representing the reference sequence. The software provides command line utilities, or it can be used as a Python library. 

# Prerequisites
`pycofitness` is implemented mainly in Python with the pseudolikelihood maximization parameter inference part implemented using C++ backend to enable parallelization. It requires: 
* Python 3
* C++ compiler that supports C++11, e.g., the GNU compiler collection (GCC)
* Optionally, OpenMP for multithreading support

# Installing
To install the current version of `pycofitness` from PyPI, run on the command line
```bash
$ pip install pycofitness
```

# Using `pycofitness` as a Python library
After installation, `pycofitness` can be imported into other Python source codes and used. For example,  

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
We can also pass parameters such as the number of iterations, the number of threads for parallel execution, and so on, as 
keyword arguments to the `PointMutation` constructor:
```python 
point_mutation = PointMutation(msa_file, biomolecule,
    max_iterations = 1000,
    num_threads = 4,
    seqid = 0.9,
    lambda_J = 10.0,
    lambda_h = 5.0,
    verbose = True
)
```
where `max_iterations` is the number of maximum iterations for gradient descent, `num_threads` is the number of 
threads for parallel execution, `seqid` is the sequence identity threshold value (if sequences have similarity more than this value, they are regarded as the same), `lambda_J` and `lambda_h` are penalizing constants for L2 regularization. If `verbose` is set to `True,` logging information is enabled.

# Running `pycofitness` from the command line
When `pycofitness` is installed, it provides a command `pycofitness` that can be executed from the command line.
For example:
```bash
$ pycofitness <biomolecule> <msa_file> --num_threads=4 --verbose
``` 
where `<biomolecule>` refers to values "protein" or "RNA" (case insensitive) and `<msa_file>` is the MSA file in the fasta format. 
The optional arguments that pycofitness input command has are listed here:

`-h, --help `: show the help message

`--seqid SEQID`: Allow the choice of the cut-off value of sequences similarity above which they are lumped together in the calculation of frequencies.

`--verbose`
`--lambda_J`
`--lambda_h`

-h, --help            show this help message and exit
  --seqid SEQID         Cut-off value of sequences similarity above which they
                        are lumped together.
  --lambda_h LAMBDA_H   Value of fields penalizing constant for L2
                        regularization of fields.
  --lambda_J LAMBDA_J   Value of couplings penalizing constant for L2
                        regularization of couplings.
  --max_iterations MAX_ITERATIONS
                        Number of iterations for gradient decent for negative
                        pseudolikelihood minimization.
  --num_threads NUM_THREADS
                        Number of threads from plmDCA computation
  --verbose             Show logging information on the terminal.
  --output_dir OUTPUT_DIR
                        Directory path to which output results are written. If
                        the directory is not existing, it will be created. If
                        this path is not provided, an output directory is
                        created using the base name of the MSA file, with
                        "output_" prefix added to it.





Information about command line options can be obtained using: 
```bash
$ pycofitness --help
```

# Preprocessing input MSAs

The input data for pycofitness is an MSA file in FASTA format. We chose not to include tools to generate MSA in pycofitness because, in this way, users can employ
their favorite MSA methods. There is only one compulsory adjustment within the input MSA: the sequence
to mutate (the reference sequence) has to be positioned as the first one in the MSA. Also, pycofitness removes MSA columns containing non-standard residues and gaps at the reference sequence.
