---
title: "Energy Tabulator Calculations"
keywords: energy
tags: [energy-tabulator]
sidebar: home_sidebar
permalink: calculations_energy-tabulator.html
summary: This page gives an overview of the energies required in the theory and how to run the EnergyTabulator code.
---

## Energies needed

There are four energies calculated in the [`EnergyTabulator`][code_energy-tabulator] code: delta function, zeroth-order, first-order, and plotting. 

The delta-function energy must represent the total energy needed to be dissipated through phonons to satisfy conservation of energy. Because we consider nonradiative transitions, this energy must be the electronic energy plus the relaxation energy. The zeroth-order term requires an electronic-only total energy difference, while the first-order term uses an electronic eigenvalue difference. We chose to plot the cross section as a function of the electronic-only total energy difference (from the zeroth-order term) in eV. 

## Inputs
The input file looks like 

```fortran
&inputParams
  ! System used for eigenvalues
  exportDirEigs = 'path-to-system-export-for-eigenvalues'

  ! Systems used for total energy differences
  exportDirInitInit = 'path-to-relaxed-initial-charge-state-export'
  exportDirFinalInit = 'path-to-final-charge-state-initial-positions-export'
  exportDirFinalFinal = 'path-to-relaxed-final-charge-state-export'

  ! Physical problem details
  captured = .true. or .false.      ! If the carrier is captured
  elecCarrier = .true. or .false.   ! If the carrier is an electron (vs a hole)

  ! Parameters for change in energy
  eCorrectTot = real				  ! size of total-energy correction in eV; default 0.0
  eCorrectEigRef = real			  ! size of correction to eig diff to ref carrier in eV; default 0.0
  refBand = integer						! band location of WZP reference carrier
  
  ! Band range for overlaps
  iBandIinit = integer						! lowest initial-state band
  iBandIfinal = integer						! highest initial-state band
  iBandFinit = integer						! lowest final-state band
  iBandFfinal = integer						! highest final-state band
  
  ! Output info
  energyTableDir = 'path-to-store-energy-tables' 			! default './'
/
```

Some things to keep in mind:
* The eigenvalues should come from a ground-state calculation because those are the most accurate. 
* The eigenvalue system should reflect the max number of k-points between the zeroth- and first-order terms.
* The total energies should ideally come from HSE calculations. 
* If PBE energies must be used, there is the ability to add corrections to the total energies and the eigenvalue of the reference state (where the WZP carrier is placed).
* The energies are always in terms of the actual electron that is transitioning (not the hole if that is the carrier).

## Outputs

An `energyTable.isp.ik` file is output for each spin and k-point in the system. The header includes the band limits for the initial and final states and the total energies used in the delta-function and zeroth-order energies.

The first two columns of the rest of the data are the indices of the initial and final bands. The delta-function, zeroth-order, first-order, and plotting energies are then all listed in order. All energies are in Hartree besides the plotting energy that is in eV.

{% include links.html %}