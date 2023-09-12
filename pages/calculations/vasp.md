---
title: "VASP"
keywords: vasp
tags: [vasp]
sidebar: home_sidebar
permalink: vasp.html
summary: This page gives a summary of the VASP calculations needed for a capture calculation, including parameter specifications and tips.
---

The following VASP calculations are needed for a capture calculation:
* relaxed perfect crystal (plus SCF after)
* relaxed final charge state (plus SCF after)
* relaxed initial charge state (plus SCF after)
* SCF for final charge state in initial positions

The displacement between the initial- and final-charge-state relaxed positions gives the $\Delta q_j$ that determines how the electronic energy is split into the different phonon modes. It is recommended to extract the energies (using [`EnergyTabulator`][calculations_energy-tabulator]) from the ground-state system in the initial relaxed positions. The wave functions are only needed from the perfect crystal and final charge state/initial positions. 

The wave functions and energies are extracted using [`Export`][calculations_export]. For any calculations where only the final positions are needed, it is not required to do an SCF calculation after or to run [`Export`][calculations_export].

## Recommended order

To run the VASP calculations, you must set up your perfect-crystal and defect-crystal supercells. The recommended order for this is
* volume relax the unit cell
* scale up to supercell
* relax inner degrees of freedom (DOF) for pristine supercell (plus SCF after)
* introduce your defect in the middle of the supercell
* relax inner DOF for defect supercell in ground state
* rest of VASP defect calculations to get all wave functions and positions needed

It is recommended to start with the ground state of the defect supercell because it is easiest to converge. The `CONTCAR`, `WAVECAR`, and `CHGCAR` from the ground-state relaxation can be used as a starting point for the rest of the defect calculations. 

If you are using more than 1 k-point with a large supercell, it is recommended to 
* relax at gamma (if relaxation needed)
* SCF at gamma
* NSCF at full k-point grid (`ICHARG=11`)

For HSE calculations, start from a preconverged PBE calculation for speed. The [`Export`][calculations_export] code is not currently implemented for the gamma-only version of VASP, except for exporting the energies only, so keep that in mind when running HSE calculations especially. If you do use the gamma-only version to get accurate HSE energies, make sure that your preconverged PBE calculation was run using the gamma-only version as well, otherwise you will get an error saying `plane wave coefficients changed`.

## Parameters

### [`EDIFF`](https://www.vasp.at/wiki/index.php/EDIFF){:target="_blank"}

For the relaxations, it is okay to use a looser `EDIFF` like `1e-5`, but we use the wave functions from the SCF calculations, so the convergence criteria must be tighter for those. We have used `1e-8` in the past. To ensure that your convergence is tight enough, run the [`TME`][calculations_tme] code with the same system input for both `PC` and `SD` and a small set of bands used for the initial state and final state. Check the `|<f|i>|^2` column in the `allElecOverlap.*.*` file to ensure that the bands are orthonormal. For the same initial and final band, the overlap should be close to 1. For different initial and final bands, the overlap should be close to zero (`~1e-10` is usually good).

### [`NBANDS`](https://www.vasp.at/wiki/index.php/NBANDS){:target="_blank"}

From the VASP wiki:
> In the electronic minimization calculations, empty states do not contribute to the total energy, however, empty states are required to achieve a better convergence. 

The default setting in VASP already includes extra bands above the filled states to get better convergence; however, we need accurate calculations of both occupied and unoccupied states when considering electron capture. For electron capture, ensure that `NBANDS` is set so that there are plenty empty bands _above the empty states that you are using for your initial electron states._ The default value from VASP is `NELECT/2+NIONS/2`. You should at least increase above the default by the number of unoccupied states that you need.

### [`ISYM`](https://www.vasp.at/wiki/index.php/ISYM){:target="_blank"}

Our code is not currently implemented with symmetrization, so symmetry must be completely turned off with `ISYM=-1`.

### [`ISPIN`](https://www.vasp.at/wiki/index.php/ISPIN){:target="_blank"}

Make sure that the proper choice of `ISPIN` is used for each of your calculations. The [`Export`][calculations_export] and [`TME`][calculations_tme] codes handle spin polarization automatically, so it is not required to use the same setting for every system. 


{% include links.html %}