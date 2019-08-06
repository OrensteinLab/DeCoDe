# DeCoDe: degenerate codon design for complete protein-coding DNA libraries

## Requirements

DeCoDe requires Python 3.7 or higher.

To run DeCoDe, you will need a local installation of [Gurobi](http://www.gurobi.com/downloads/download-center) with an appropriate lisence (academic licenses are provided for free direct from Gurobi). You will also need to install [`gurobipy`](https://www.gurobi.com/documentation/8.1/quickstart_mac/the_gurobi_python_interfac.html) and make sure that you are able to import it from within your local environment.

Finally, you will need the custom version of CVXPY included in this repository. Due to a [known issue](https://github.com/cvxgrp/cvxpy/issues/735), CVXPY does not currently directly support the use of time limits with the Gurobi solver. The custom CVXPY provided here patches this issue and can be installed by running `pip install .` from inside the subdirectory.

## Setup

1. Follow the directions to install Gurobi outlined [here](https://www.gurobi.com/documentation/8.1/quickstart_mac/the_gurobi_python_interfac.html). If appropriate, [request an academic license](https://www.gurobi.com/documentation/8.1/quickstart_linux/obtaining_a_gurobi_license.html).
2. Clone this directory: `git clone https://github.com/OrensteinLab/DeCoDe.git`
3. Change into the cloned directory: `cd DeCoDe`
4. If you would like to run DeCoDe under a time constraint, install the local version of `CVXPY` found in this directory:
```bash
cd cvxpy
pip install .
```
5. Change back into the parent directory: `cd ..`
6. Run `pip install -r requirements.txt` to install all remaining requirements.

## Testing your setup

To test your setup, you can run the following command which will optimize the small example library found in Figure 1B of the paper describing DeCoDe:

```
python decode.py --limit 6 \
    --sublib 1 \
    --threads 10 \
    examples/fig1/6o_test-3.aln \
    test.json
```

Proper output will look something like this:

```
Reading MSA file...
Removing fixed positions...

Library size limit:		6
Sublibrary count limit:		1
Number of targets:		6
Number of variable positions:	5

Using exact library size.

Academic license - for non-commercial use only
Parameter OutputFlag unchanged
   Value: 1  Min: 0  Max: 1  Default: 1
Changed value of parameter QCPDual to 1
   Prev: 0  Min: 0  Max: 1  Default: 0
Changed value of parameter Threads to 10
   Prev: 0  Min: 0  Max: 1024  Default: 0
Optimize a model with 30 rows, 4217 columns and 25398 nonzeros
Variable types: 0 continuous, 4217 integer (4217 binary)
Coefficient statistics:
  Matrix range     [7e-01, 5e+00]
  Objective range  [1e+00, 1e+00]
  Bounds range     [1e+00, 1e+00]
  RHS range        [1e+00, 5e+00]
Found heuristic solution: objective 0.0000000
Presolve removed 18 rows and 4194 columns
Presolve time: 0.01s
Presolved: 12 rows, 23 columns, 86 nonzeros
Variable types: 0 continuous, 23 integer (23 binary)

Root relaxation: objective -4.600000e+00, 21 iterations, 0.00 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0   -4.60000    0    3    0.00000   -4.60000      -     -    0s
H    0     0                      -3.0000000   -4.60000  53.3%     -    0s
     0     0     cutoff    0        -3.00000   -3.00000  0.00%     -    0s

Cutting planes:
  Gomory: 2
  StrongCG: 1
  GUB cover: 3

Explored 1 nodes (32 simplex iterations) in 0.01 seconds
Thread count was 10 (of 40 available processors)

Solution count 2: -3 0

Optimal solution found (tolerance 1.00e-04)
Best objective -3.000000000000e+00, best bound -3.000000000000e+00, gap 0.0000%

Number of covered targets:	3
Total library size:		6
Probability on target:		0.50000

Writing output...
```

# Usage

DeCoDe requires input sequences to be pre-aligned in the ClustalW format. To generate the the alignment file (`.aln`), you can use the Clustal Omega server found [here](https://www.ebi.ac.uk/Tools/msa/clustalo/) or through any other program that supports output to the ClustalW file format.

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

A full help page is available through the `python decode.py --help` command, the output of which is reproduced below:

```
Usage: decode.py [OPTIONS] ALIGNMENT_FILE OUTPUT_FILE

  Optimize a degenerate codon library given a set of target sequences, a
  total size limit, and a sublibrary count limit.

Options:
  --limit INTEGER       Total library size limit.  [required]
  --sublib INTEGER      Total sublibrary limit  [required]
  --bins INTEGER        Specify the number of bins for approximation of multi-
                        sublibrary optimizations.  [default: 100]
  --time-limit INTEGER  Time limit in seconds for the ILP solver.  [default:
                        0]
  --threads INTEGER     Time limit in seconds for the ILP solver.  [default:
                        0]
  -q, --quiet           Run quietly.  [default: False]
  --help                Show this message and exit.
```

# Interpreting output:

Below is the output file (`test.json`) from the following command:

```
python decode.py --limit 6 \
    --sublib 1 \
    --threads 10 \
    examples/fig1/6o_test-3.aln \
    test.json
```

Output:

```
{
    "solution_optimal": true,
    "fixed_positions": {
        "0": "D",
        "1": "G",
        "2": "D",
        "4": "N",
        "6": "H",
        "8": "F",
        "10": "V",
        "12": "G",
        "13": "E"
    },
    "variable_positions": [
        3,
        5,
        7,
        9,
        11
    ],
    "sequences": [
        "VGKSS",
        "VVKSS",
        "VARSS",
        "AGKSS",
        "AV-TT",
        "AA-TT"
    ],
    "coverage": [
        1.0,
        1.0,
        0.0,
        1.0,
        0.0,
        0.0
    ],
    "n_var_pos": 5,
    "n_covered": 3,
    "total_lib_size": 6,
    "on_target_p": 0.5,
    "parsed_lib": [
        [
            "D",
            "G",
            "D",
            {
                "A/V": [
                    "GYA",
                    "GYC",
                    "GYG",
                    "GYT"
                ]
            },
            "N",
            {
                "A/G/V": [
                    "GBA",
                    "GBC",
                    "GBG",
                    "GBT"
                ]
            },
            "H",
            "K",
            "F",
            "S",
            "V",
            "S",
            "G",
            "E"
        ]
    ],
    "construct_time": 0.16756463050842285,
    "solve_time": 0.011950969696044922,
    "total_time": 0.17951560020446777
}
```

Here, we define all of the above entries in the JSON file:

- `solution_optimal` - Whether the found solution is optimal.
- `fixed_positions` - All positions and residue identities fixed across the set of target sequences.
- `variable_positions` - All positions that vary across the target set.
- `sequences` - All of the target sequences.
- `coverage` - Whether the corresponding target in `sequences` is covered by the produced library.
- `n_var_pos` - Total number of variable positions in the aligned sequences.
- `n_covered` - Total number of full length target sequences covered by the library (= `sum(coverage)`).
- `total_lib_size` - Total produced library size.
- `on_target_p` - Fraction of the produced library that covers a target sequence (= `n_covered / total_lib_size`)
- `parsed_lib` - A list of the individual sublibraries. Each list contains the one letter code of the amino acid at that position if the position is fixed in the designed library. If the position is variable and DeCoDe has allocated a degenerate codon for the given position, the output key will include a `/`-separated list of the covered amino acids and a list of all equivalent degenerate codons from which the user can choose a codon to employ in the finished library.
- `construct_time` - The total time (in seconds) for CVXPY to construct the problem and hand it off to the Gurobi solver.
- `solve_time` - Total time for Gurobi to solve the design problem.
- `total_time` - Total time = `construct_time` + `solve_time`.

## Running library designs in parallel

We can use [GNU `parallel`](https://www.gnu.org/software/parallel/) to design several libraries in parallel. These examples will generate all of the libraries for the demonstrations present in the manuscript.

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
parallel --will-cite --lb \
    "/usr/bin/time -v --output results/sl_comparison/sl/gfp_239_{1}_{2}.time \
    python run_SwiftLib.py --limit {1} --sublib {2}" \
    ::: 100000 1000000 10000000 100000000 1000000000 ::: 1 2
```
