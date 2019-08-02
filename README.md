# DeCoDe: degenerate codon design for complete protein-coding DNA libraries

## Requirements

DeCoDe requires Python 3.7 or higher.

To run DeCoDe, you will need a local installation of [Gurobi](http://www.gurobi.com/downloads/download-center) with an appropriate lisence (academic licenses are provided for free direct from Gurobi). You will also need to install [`gurobipy`](https://www.gurobi.com/documentation/8.1/quickstart_mac/the_gurobi_python_interfac.html) and make sure that you are able to import it from within your local environment.

Finally, you will need the custom version of CVXPY included in this repository. Due to a [known issue](https://github.com/cvxgrp/cvxpy/issues/735), CVXPY does not currently directly support the use of time limits with the Gurobi solver. The custom CVXPY provided here patches this issue and can be installed by running `pip install .` from inside the subdirectory.

## Setup

1. Follow the directions to install Gurobi outlined [here](). If appropriate, request an academic license.
2. Clone this directory: `git clone https://github.com/OrensteinLab/DeCoDe.git`
3. If you would like to run DeCoDe under a time constraint, install the local version of `CVXPY` found in this directory:
```bash
pip install DeCoDe/cvxpy
```

# Usage

Change into the DeCoDe directory and run the algorithm using a variant of the command below, replacing placeholders with your desired files where appropriate:

```
python decode.py \
	--limit <library limit> \
	--sublib <number of sublibraries> \
	--bins <number of bins, if applicable> \
	--time-limit <time limit, in seconds> \
	--threads <number of threads for ILP solver parallelization> \
	<input file>.aln \
	<output file>.json
```

## Running library designs in parallel

We can use [GNU `parallel`](https://www.gnu.org/software/parallel/) to design several libraries in parallel.

```bash
parallel --will-cite --lb -j 3 \
    "/usr/bin/time -v --output results/multi_sublib/gfp_exclude_long_10000000_{1}.time \
    python decode.py --limit 10000000 --sublib {1} --bins 100 --time-limit 172800 \
    --threads 12 examples/gfp/gfp_exclude_long.aln \
    results/multi_sublib/gfp_exclude_long_10000000_{1}.json | \
    tee results/multi_sublib/gfp_exclude_long_10000000_{1}.log" ::: 1 2 3 4 8 12
```

```bash
parallel --will-cite --lb -j 3 \
    "/usr/bin/time -v --output results/sl_comparison/ilp/gfp_239_{1}_{2}.time \
    python decode.py --limit {1} --sublib {2} --bins 100 --time-limit 172800 \
    --threads 12 examples/gfp/gfp_239.aln results/sl_comparison/ilp/gfp_239_{1}_{2}.json | \
    tee results/sl_comparison/ilp/gfp_239_{1}_{2}.log" \
    ::: 100000 1000000 10000000 100000000 1000000000 ::: 1 2
```

```bash
parallel --will-cite --lb "/usr/bin/time -v --output results/sl_comparison/sl/gfp_239_{1}_{2}.time python run_SwiftLib.py --limit {1} --sublib {2}" ::: 100000 1000000 10000000 100000000 1000000000 ::: 1 2
```
