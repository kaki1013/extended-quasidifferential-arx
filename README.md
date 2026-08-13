# Extended Quasidifferential Analysis of SPECK and CHAM

This repository contains the research code accompanying the manuscript
**"Revisiting Differential Characteristics of SPECK and CHAM via Extended
Quasidifferential Framework,"** submitted to *IEICE Transactions on
Fundamentals of Electronics, Communications and Computer Sciences*.

The manuscript is currently under submission. This repository does not imply
acceptance or publication by IEICE.

## Overview

The code adapts the extended quasidifferential framework to the ARX block
ciphers SPECK and CHAM. It builds SMT models that enumerate quasidifferential
trails compatible with a fixed differential characteristic while incorporating
the cipher's key schedule. The signed trail contributions can then be summed
to study the expected differential probability (EDP) relative to the usual
independence-based estimate.

The repository also contains the key-schedule variants used in the manuscript
to compare nonlinear and linear schedules while keeping the SPECK round
function fixed.

## Repository layout

| Path | Description |
| --- | --- |
| `speck_original/` | SPECK round function with the original nonlinear SPECK key schedule. |
| `speck_linear_key_sch_identity/` | SPECK round function with the simple cyclic linear schedule, where round `i` uses `K[i mod m]`. |
| `speck_linear_key_sch_cham-like/` | SPECK round function with the linear key schedule derived from CHAM. |
| `cham/` | CHAM analysis with its original linear key schedule. |
| `zero_verification/` | Exhaustive C++ checks used for the zero-correlation verification discussed in the manuscript. |

Within each Python analysis directory:

- `common.py` constructs the SMT model, enumerates solutions, and computes
  correlation signs.
- `utils.py` implements bit-vector operations and the modular-addition QDTM
  constraints.
- `speck.py` or `cham.py` analyzes one selected characteristic.
- `speck_all_data.py` or `cham_all_data.py` processes all matching
  characteristic files in the local `data/` directory.

## Experimental environment

The experiments reported in the manuscript used:

- Ubuntu 24.04.1 LTS
- Python 3.10.20
- Boolector 3.2.4 through PyBoolector 3.2.4
- GCC with C++17 support for the exhaustive verification programs

PyBoolector must be importable as `pyboolector`. For example, in a suitable
Python environment:

```bash
python -m pip install pyboolector==3.2.4
```

Exact package availability may depend on the operating system and Python
distribution. Boolector can alternatively be built with its Python bindings
enabled.

## Characteristic input format

Characteristic files are read from a `data/` directory below the selected
implementation directory. Blank lines and lines beginning with `#` are
ignored.

SPECK files use one pair of hexadecimal word differences per line:

```text
left_difference,right_difference
```

CHAM files use four hexadecimal word differences per line:

```text
x0_difference,x1_difference,x2_difference,x3_difference
```

Each file contains the state differences at the round boundaries, so an
`r`-round characteristic contains `r + 1` lines. The expected filenames are:

```text
speck_<block>_<rounds>r_<weight>_<index>.txt
cham_<block>_<rounds>r_<weight>_<index>.txt
```

The characteristic datasets are research inputs and should be placed in the
corresponding local `data/` directory before running the scripts.

## Running a single characteristic

First edit the parameter block near the top of `speck.py` or `cham.py` to
select the block size, round count, characteristic weight, index, key size or
number of key words, and enumeration cutoff. Then run the script from its own
directory so that relative paths resolve correctly.

For example:

```bash
cd speck_original
python speck.py
```

or:

```bash
cd cham
python cham.py
```

The SPECK script prints the characteristic, enumerated mask trails, their
signs and relative weights, and the cumulative normalized sum. The CHAM script
writes analogous summaries to a result file when `save_file = True`.

## Running all available characteristics

From the desired implementation directory, run:

```bash
python speck_all_data.py
```

or:

```bash
python cham_all_data.py
```

These scripts scan the local `data/` directory and evaluate every filename
matching the format above for all supported key lengths.

## Enumeration cutoff and reported values

For a quasidifferential trail, the code uses `w` for the correlation-weight
loss relative to the all-zero-mask trail. Its normalized signed contribution
is

```text
sign * 2^(-w).
```

`max_weight_loss` is an exclusive upper bound: the scripts enumerate
`0 <= w < max_weight_loss`. In the notation of the manuscript, choosing
`max_weight_loss = t` retains trails whose relative-correlation magnitude is
at least `2^(-(t-1))`; adjust the bound when reproducing a particular table or
figure.

The accumulated value printed by the current scripts is the normalized signed
sum, corresponding to the relative EDP (REDP). Multiplying it by the product
of the one-round differential probabilities gives the truncated EDP estimate.
Because trails below the cutoff are omitted, the result is an estimate rather
than an exact EDP unless the enumeration is complete.

## Exhaustive zero-correlation checks

The programs in `zero_verification/` exhaustively test the two modular-addition
transitions used in the manuscript's zero-correlation argument.

```bash
cd zero_verification
make
./41_upper
./41_lower
```

Each program examines all `2^32` pairs of 16-bit inputs. Consequently, runtime
can be very long without additional optimization or parallelization.

## Reproducibility notes

- Run each script from the directory containing it; input and output paths are
  relative to the current working directory.
- Record the chosen parameter block and `max_weight_loss` with every result.
- Enumeration time depends strongly on the characteristic, block size, key
  schedule, and cutoff.
- The code evaluates supplied characteristics; it is not a general-purpose
  differential-characteristic search tool.

## Acknowledgments

Portions of `utils.py` and `common.py` are used from or adapted from Tim
Beyne's [`quasidifferential-trails`](https://github.com/TimBeyne/quasidifferential-trails)
repository. In particular, that repository provided the basis for parts of
the bit-vector utilities and SMT modeling of quasidifferential trails. We
gratefully acknowledge Tim Beyne's original implementation.

The corresponding quasidifferential framework is described in:

> T. Beyne and V. Rijmen, "Differential Cryptanalysis in the Fixed-Key Model,"
> in *Advances in Cryptology - CRYPTO 2022*, LNCS 13509, pp. 687-716, 2022.
> DOI: [10.1007/978-3-031-15982-4_23](https://doi.org/10.1007/978-3-031-15982-4_23).

## Citation

As the manuscript is still under submission, final bibliographic metadata is
not yet available. Until publication information is assigned, please cite it
as an unpublished manuscript:

```bibtex
@unpublished{lee2026revisiting,
  author = {Myungkyu Lee and Yunjae Hwang and Hanbeom Shin and Insung Kim and
            Sunyeop Kim and Byoungjin Seok and Dongjae Lee and Deukjo Hong and
            Jaechul Sung and Seokhie Hong},
  title  = {Revisiting Differential Characteristics of SPECK and CHAM via
            Extended Quasidifferential Framework},
  note   = {Manuscript submitted to IEICE Transactions on Fundamentals},
  year   = {2026}
}
```

## License

This project is available under the [MIT License](LICENSE). Copyright (c) 2026
Myungkyu Lee.
