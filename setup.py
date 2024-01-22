from setuptools import setup, find_packages, Extension
#from setuptools.command.build_ext import build_ext 

with open("README.md") as fh:
    long_description = fh.read()

requirements = [
    'llvmlite>0.30.0',
    'numpy',
    'biopython>=1.74',
]

plmdca_compile_args = ["-fopenmp", "-std=c++11", "-O3"]  
plmdca_link_args = ["-fopenmp", "-O3"] 


plmdca_ext = Extension(
    'pycofitness.plmdca._plmdcaBackend',
    [   
        'pycofitness/plmdca/lbfgs/lib/lbfgs.cpp',
        'pycofitness/plmdca/plmdca_numerics.cpp',
        'pycofitness/plmdca/plmdcaBackend.cpp', 
    ],
    include_dirs=[
        'pycofitness/plmdca/include/',
        'pycofitness/plmdca/lbfgs/include/',
    ],
    extra_compile_args = plmdca_compile_args,
    extra_link_args = plmdca_link_args,
    language = "c++",  
)

setup(
    name="pycofitness",
    version="1.3",
    author="Fabrizio Pucci, Mehari B. Zerihun",
    author_email="Fabrizio.Pucci@ulb.be, mbzerihun@gmail.com",
    python_requires=">=3.6",
    description="<i>In silico</i> mutagenesis of protein and RNA sequences using coevolution",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KIT-MBS/pycofitness",
    download_url="https://pypi.org/project/pycofitness/",
    packages=find_packages(
        exclude=["*.tests","*.tests.*","tests.*", "tests",
            "examples", "*.examples", "examples.*", "*.examples.*",
        ],
    ),
    ext_modules = [plmdca_ext],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "Programming Language :: C",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 4 - Beta",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires= requirements,
    tests_require = requirements,
    entry_points={
        "console_scripts":[

            "pycofitness=pycofitness.main:run_mutation",
        ],
    },
    test_suite="tests",
)
