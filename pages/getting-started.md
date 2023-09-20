---
title: "Getting started with the Defect Cross Sections code"
keywords: setup start
tags: [getting_started]
sidebar: home_sidebar
permalink: getting-started.html
summary: These brief instructions will help you get started quickly with the code. The other topics in the code go more in-depth into theory, the code itself, and examples.
---

## Setup

First, you will need to have VASP set up on the machine you are working on. If it is already compiled on the machine, it is best to use that executable because it will be well-optimized for the specific machine structure. If not, work with the help desk to set up the required `makefile.include` to get VASP compiled. 

{% include note.html content="The [`Export`][calculations_export] currently only works with the standard version of VASP or to export only the energies from results from the gamma-only version of VASP. Assumptions made in the code are not compatible with the noncollinear version." %}

To setup the defect cross sections code, clone the repo, then edit the parameters at the top of the main `Makefile`. The QE-based versions of `Export` is currently out of date, so don't worry about updating those values. Make sure to update the compilers `f90` and `mpif90` to the intel wrappers specific to your system.

{% include note.html content="The `Export` and `TME` codes utilize binary I/O for parallelization. This feature has only been tested with an Intel compiler, so it may break with another compiler." %}

To compile all of the code needed to run with VASP, run `make all_VASP`. Running `make all` will currently fail because the QE-based codes are still intertwined with the QE code itself. We have not had a need to implement the QE version right now because we ran into issues with doing HSE calculations with d-electrons in QE, so it has not been a focus. The feature can always be implemented in the future if needed.

## Running the codes

The [home page][index] contains an overview of how the codes piece together, and the [theory][theory] page contains a brief overview of the theory and how it relates to the codes. There are many choices to make in the VASP codes, and care must be taken to pass the right inputs into the codes in the correct order. It is _highly recommended_ to review those pages before running the codes.

Once you have reviewed the overview pages to get a sense of how the codes work together, check out the individual documentation pages for detailed information on how to run each of the codes:
* [`VASP`][vasp]
* [`Export`][calculations_export]
* [`EnergyTabulator`][calculations_energy-tabulator]
* [`Phonons`][phonons]
* [`PhononPP`][calculations_phononpp]
* [`TME`][calculations_tme]
* [`LSF`][calculations_lsf]

## Tutorial

As you are learning, it is helpful to run an example calculation. Check out the [tutorial][tutorial] page on running a capture calculation on a triply-hydrogenated Si vacancy.

{% include links.html %}
