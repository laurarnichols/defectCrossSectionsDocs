---
title: Introduction
sidebar: home_sidebar
permalink: index.html
keywords: homepage introduction
tags: [introduction]
---

## Overview

This site provides the documentation for the [Defect Cross Sections project](https://github.com/laurarnichols/defectCrossSections){:target="_blank"}. This code currently calculates the capture cross section for a defect in a semiconductor based on the [theory][theory] developed in the Barmparis and GaN papers. There are codes to integrate with both VASP and QuantumEspresso, but __only the VASP version is currently up-to-date__.

The [flow chart for the full calculation](#complete-process) may be overwhelming, but _stick with me_. It is important to get an overall understanding of how the code pieces together because there are so many moving parts. But before showing you the full flow chart, let's start with a simple one:

{% include image.html file="CaptureFlowChartSimple.png" width="40%" alt="Simple capture flow chart" %}

<br/> 

The basic process for running a capture calculation is to get all of the VASP calculations that you need done, run `Export` to translate the information into a form easily read by the rest of the code, calculate phonons, post-process the files to gather various energies needed and information about the phonons, then run `TME` and `LSF` for the zeroth- and first-order terms. 

## VASP and Export Calculations 

The following VASP calculations are needed for a capture calculation:
* relaxed perfect crystal (plus SCF after)
* relaxed final charge state (plus SCF after)
* relaxed initial charge state (plus SCF after)
* final charge state in initial positions

See the [VASP][vasp] page to find details on how the calculations should be set up and the suggested order of running the calculations.

The results must then be post-processed by the `Export` code to be used as input to the other programs in the suite. Once the VASP results are exported, the raw files are not needed again. Some of the items exported are:
* wave functions
* projectors (for PAW)
* projections
* plane-wave grid
* pseudopotential

See the [`Export`][calculations_export] page for more details on how to run the `Export` code and what outputs to expect.

## Energy tabulation

There are difference energies used in different parts of the theory. Instead of calculating the different energies in different places, the `EnergyTabulator` code is used to gather the energy information from the exported data from the various systems and tabulate all of the needed energies for use in the `TME` and `LSF` codes. See the [`EnergyTabulator`][calculations_energy-tabulator] page for more details.

## Phonons and phonon post-processing

The phonon eigenvectors and frequencies are needed to determine how the electronic-transition energy is dissipated into each of the phonon modes. It is assumed that eigenvectors and frequencies of the phonons in the initial and final relaxed positions are the same. Consequently, it does not matter which equilibrium positions are used for the phonon calculations; however, it is recommended to use the ground electronic state of the system for faster convergence of the phonon calculations. 

The phonon post-processing code (`PhononPP`) is then used to generate positions shifted along each eigenvector (needed for the first-order matrix element) and the Huang-Rhys factor (needed for distributing the energy across the modes in the `LSF` code). `PhononPP` is currently set up to read the required information in the form of the `mesh.yaml` output by `Phonopy`. 

Check out the [Phonons][phonons] and [`PhononPP`][calculations_phononpp] pages for guidance on running a phonon calculation with `Phonopy` and how to use `PhononPP`.

## Zeroth-order

{% include note.html content="The zeroth-order term is only valid for nonequilibrium capture. See the [theory][theory] page for more details." %}

The zeroth-order transition rate is straightforward to calculate: run the `TME` and `LSF` codes with `order=0`. The [`TME`][calculations_tme] code calculates the transition matrix element, then the [`LSF`][calculations_lsf] code does the time-domain integral to get the final transition rate. 

## First-order

The first-order matrix elements are more complex because a matrix element $M_j$ must be calculated for each phonon mode $j$. The [`PhononPP`][calculations_phononpp] code generates `POSCAR` files with a small displacement along each of the phonon eigenvectors. The wave functions are needed for each of the shifted `POSCAR`s, so a VASP SCF calculation must be run, then the data must be post-processed using `Export`. The `TME` code is then run with the wave functions from the shifted postitions and the wave functions from the unshifted positions.

{% include note.html content="We have found that the charge state must be the same for both cases, otherwise $M_j$ will not be independent of the displacement $\delta q_j$. The charge state used will not significantly affect the wave functions, but it must be the same for both sets of inputs. The positions must be the initial-charge-state equilibrium positions (and small displacements from them)." %}

Once the first-order matrix elements are calculated, run the `LSF` code to get the first-order transition rate. 

## Complete process

Altogether, the complete process for a capture calculation is:

{% include image.html file="CaptureFlowChart.png" width="70%" alt="Capture flow chart" %}

[Jump back to the top](#overview) to get a breakdown of each piece. 

## Getting started

To get started, see [Getting Started][getting-started].

{% include links.html %}
