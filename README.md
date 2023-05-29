# About `pycofitness`
`pycofitness` is Python implemementation of *in silico* mutagenesis for multipe sequence alignments (MSA) of protein or RNA sequence families using **__pseudolikelihood maximization__**  direct couplings analysis (DCA) algorithm. The input file for `pycofitness` must be an MSA of sequences in fasta format, with the first sequence representing the reference sequence. The software provides command line utilities or it can be used as a Python library. 

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

# Using `pycofitness` as a Python Library
After installation, `pycofitness` can be imported into other Python source codes and used. For example,  

```python 
from pycofitness.mutation import PointMutation

point_mutation = PointMutation(msa_file, biomelcule)
deltas = point_mutation.delphi_epistatic()
#print the values of delta_dict by iteration over sites and bases/residues

for i in range(len(deltas)):
    for res in deltas[i]:
        print(i + 1, res, deltas[i][res])
#Note we added 1 to i to count sites starting from one as opposed to zero.
```
We can also pass parameters such as number of iterations, and number of threads for parallel execusion, and so on, as 
keyword arguments to `PointMutation` constructor:
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
where `max_iterations` is the number of maximum iterations for gradient decent, `num_threads` is the number of 
threads for parallel execution, `seqid` is sequence identity cut off value (if sequences have similarity more that this value, they are regarded as the same), `lambda_J` and `lambda_h` are penalizing constants for L2 regularization. If `verbose` is set to be `True`, logging information is enabled.

# Running `pycofitness` From Command Line
When `pycofitness` is installed, it provides a command `pycofitness` that can be executed from the command line.
For example:
```bash
$ pycofitness <biomolecule> <msa_file> --num_threads=4 --verbose
``` 
where `<biomolecule>` refers values "protein" or "RNA" (case insensitive) and `<msa_file>` is the MSA file. 

Information about command line options can be obtained using: 
```bash
$ pycofitness --help
```
