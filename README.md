# DeCoDe: degenerate codon design for complete protein-codingDNA libraries

## Requirements

DeCoDe requires Python 3.7 or higher.

To run DeCoDe, you will need a local installation of [Gurobi](http://www.gurobi.com/downloads/download-center) with an appropriate lisence (academic licenses are provided for free direct from Gurobi). You will also need to install [`gurobipy`](https://www.gurobi.com/documentation/8.1/quickstart_mac/the_gurobi_python_interfac.html) and make sure that you are able to import it from within your local environment.

Finally, you will need the custom version of CVXPY included in this repository. Due to a [known issue](https://github.com/cvxgrp/cvxpy/issues/735), CVXPY does not currently directly support the use of time limits with the Gurobi solver. The custom CVXPY provided here patches this issue and can be installed by running `pip install .` from inside the subdirectory.

## Running library designs in parallel

We can use [GNU `parallel`](https://www.gnu.org/software/parallel/) to design several libraries in parallel.

```
parallel --will-cite --lb -j 3 "/usr/bin/time -v --output results/multi_sublib/gfp_exclude_long_10000000_{1}.time python decode.py --limit 10000000 --sublib {1} --bins 100 --time-limit 172800 --threads 12 examples/gfp/gfp_exclude_long.aln results/multi_sublib/gfp_exclude_long_10000000_{1}.json | tee results/multi_sublib/gfp_exclude_long_10000000_{1}.log" ::: 1 2 3 4 8 12
```