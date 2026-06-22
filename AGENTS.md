# Repository Guidelines

## Project Structure & Module Organization

This repository contains numerical simulation examples for computational physics.
Top-level directories are organized by method and topic:

- `FDS/`: finite difference method examples.
- `FVS/`: finite volume method examples.
- `Vlasov/`: kinetic Vlasov-Poisson simulations.
- `common/`: shared C++ utilities and Python visualization helpers.

Most topic directories are split into `C++/` and `Py/` implementations; under
`Vlasov/`, this split appears inside `Elesta/` and `Self_g/`. C++ task
directories contain `main.cpp`, a local solver implementation such as
`adv1d.cpp`, `dif1d.cpp`, or `poi2d.cpp`, a `cnst.hpp` parameter file, and a
`Makefile`. Python task directories contain `main.py` batch scripts and, where
available, notebooks for Colab. C++ example figures are stored under `imgs/`.

## Build, Test, and Development Commands

Build shared C++ code first:

```sh
cd common
make
```

Build and run an individual C++ task:

```sh
cd FDS/C++/Adv
make
./a.out
```

Remove C++ build products in a task directory:

```sh
make clean
```

Run a Python version from its task directory:

```sh
cd FDS/Py/Adv
python main.py
```

There is no central test suite. Use successful compilation and representative
program runs as the primary verification.

## Coding Style & Naming Conventions

C++ code uses compact procedural examples with two-space indentation. Keep
constants in `cnst.hpp`, shared boundary/output routines in `common/`, and
task-specific numerical schemes in local `*.cpp`/`*.hpp` pairs. Preserve the
existing short numerical names, such as `adv1d`, `dif1d`, `bgs1d`, and `poi2d`.

Python examples provide `main.py` batch scripts and, in many directories,
matching notebooks. Follow the existing lowercase variable names and keep
user-set parameters near the top of the file.

## Working rules

Keep changes minimal and easy to review.

Do not change numerical algorithms, discretization schemes, boundary conditions, or physical models unless explicitly requested.

Do not change input/output file formats unless explicitly requested.

Do not introduce external dependencies without asking first.

Do not change compiler flags, OpenMP settings, MPI settings, floating-point precision, or optimization flags without explaining the reason and risk.

Prefer explaining the plan before editing when the task affects core numerical routines.

## Testing Guidelines

Before submitting C++ changes, run `make clean && make` in the affected task
directory. If shared code changes, also run `make` in `common` and rebuild at
least one representative task. For Python changes, run the affected `main.py`
and check that the expected plot appears without exceptions.

## Done when

The code builds successfully, or the remaining failure is clearly explained.

The diff is small and understandable.

Any numerical, performance, or compatibility risks are explicitly noted.

## Commit & Pull Request Guidelines

Recent commits use short imperative summaries, for example `update Makefile by
codex` or `add poi2d_gs by codex`. Keep commit messages concise and mention the
affected component. Pull requests should describe the numerical problem touched,
the commands run for verification, and any expected output changes such as new
`.dat` files or figures.
