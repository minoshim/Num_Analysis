# Implementation note: advection equation solver from the paper

## Purpose

This document describes a task for Codex: implement a new numerical solver for the advection equation based on the algorithm proposed in the accompanying paper.

The implementation should be added to the existing `Num_Analysis` repository as a small, reviewable extension. The goal is not to redesign the repository, but to add one new solver example consistent with the existing style of the repository.

## Repository context

This repository contains numerical analysis examples for computational physics, including finite-difference, finite-volume, Poisson, and Vlasov solvers.

Relevant directories include:

* `FDS/`: finite-difference method examples
* `FVS/`: finite-volume method examples
* `Vlasov/`: kinetic / Vlasov simulation examples
* `common/`: common C++ utility code

The repository contains both C++ and Python examples. For this task, focus first on the C++ implementation unless otherwise instructed.

## Reference material

The paper PDF is placed in:

```text
references/
```

Before implementing, read the paper and identify the specific algorithm for the advection equation solver.

If a text-extracted version of the paper exists, prefer reading it first:

```text
references/*.txt
```

If there is ambiguity between the PDF and extracted text, report the ambiguity before implementing.

## Target task

Implement the advection equation solver proposed in the paper.

The target equation is the scalar advection equation, for example:

```math
\frac{\partial f}{\partial t} + v \frac{\partial f}{\partial x} = 0
```

or the corresponding form used in the paper.

The implementation should be added as a new example, rather than replacing existing solvers.

## Initial implementation scope

Implement the simplest meaningful version first:

* one-dimensional scalar advection equation
* uniform grid
* periodic boundary condition, unless the paper requires another condition
* constant advection velocity
* one or two standard initial conditions, such as:

  * Gaussian pulse
  * square wave
  * sine wave
* output format consistent with existing examples in the repository

Do not implement multi-dimensional extensions, nonlinear systems, AMR, complicated boundary conditions, or performance optimizations unless explicitly requested later.

## Required pre-implementation analysis

Before editing any source files, do the following:

1. Inspect the existing repository structure.
2. Identify the closest existing advection solver example.
3. Read the relevant part of the paper.
4. Summarize the algorithm in implementation-oriented form.
5. Identify which equations, reconstruction formulas, flux formulas, or time-integration steps are required.
6. Propose the exact files to add or modify.
7. Propose how to build and run the new example.
8. Propose a minimal verification plan.

Do not edit files until this analysis is complete.

## Constraints

* Do not migrate the project to CMake.
* Use the existing Makefile-based style.
* Keep changes minimal and easy to review.
* Do not rewrite unrelated examples.
* Do not change existing solver behavior.
* Do not change common utilities unless necessary.
* If common utilities must be changed, explain why first.
* Do not introduce external dependencies.
* Do not change file formats used by existing examples unless explicitly requested.
* Preserve the educational style of the repository.

## Numerical implementation requirements

When implementing the solver:

* keep the code readable rather than overly optimized;
* use variable names consistent with existing examples where possible;
* clearly separate:

  * grid setup,
  * initial condition,
  * numerical flux or update formula,
  * time stepping,
  * output;
* add comments for the algorithmic steps, but do not copy long passages from the paper;
* explicitly state the CFL condition used;
* avoid hidden changes to the numerical method.

## Verification requirements

After implementation, run the appropriate build command for the new example.

For C++ examples, follow the existing repository style. Typically:

```bash
make
./a.out
```

If the example depends on files in `common/`, first ensure that `common` builds correctly:

```bash
cd common
make clean
make
```

Then return to the new example directory and build it.

The verification should include at least one of the following:

1. periodic advection of a smooth profile over one period;
2. comparison between initial and final profiles after one full advection cycle;
3. qualitative comparison against an existing solver in the repository;
4. convergence check if easy to add.

If a full verification is too large for the first implementation, add a short note explaining what remains to be checked.

## Suggested output

The new example should output data in a format that can be plotted with existing Python scripts, if possible.

If a new plotting script is needed, add a minimal Python script that reads the output and plots the numerical solution.

Do not add notebook files in the first implementation unless explicitly requested.

## Done criteria

The task is done when:

* the new solver example is added;
* the code builds with the existing Makefile-based workflow;
* the example runs successfully;
* output data are generated;
* the algorithmic correspondence to the paper is summarized;
* limitations and uncertain points are listed;
* the git diff is small and reviewable.

## First Codex instruction

Start by reading this file, the repository README, the existing advection examples, and the paper in `references/`.

Then produce an implementation plan.

Do not modify files yet.
