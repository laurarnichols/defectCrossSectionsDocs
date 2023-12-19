---
title: "VASP and Export"
keywords: vasp export
tags: [vasp export]
sidebar: home_sidebar
permalink: vasp-and-export.html
summary: This page gives a summary of the VASP calculations needed for a capture calculation, including parameter specifications and tips, followed by details on how to Export those results for input into subsequent calculations.
---

## Overview 

The following systems are needed for a capture calculation:
* relaxed perfect crystal
* relaxed final-charge-state defect
* relaxed initial-charge-state defect 
* final-charge-state defect in initial positions

We need to extract the total energy from all of the defect calculations, and we need the eigenvalues from the ground state of the initial positions (ground state is most accurate). We need to extract the wave functions from the perfect crystal and the final charge state/initial positions for the zeroth-order matrix element. We also need the displacement between the initial- and final-charge-state relaxed positions to get the $\Delta q_j$ that determines how the electronic energy is split into the different phonon modes. 

{% include note.html content="It is easy enough to do the above calculations at the HSE level, but the first-order term is not currently feasible at the HSE level in large supercells (i.e., ~1500 calculations $\times$ ~7 hours per calculation for a 4x4x4 supercell on 768 processors). Both terms should be done at the same level, so the wave functions for the zeroth-order term should not be done at the HSE level at this point." %} 

The wave functions and energies are extracted using [`Export`](#export). For any calculations where only the final positions are needed, it is not required to do an SCF calculation after or to run [`Export`][#export].

### Total energies and eigenvalues

Capture involves promoting carriers between band states and defect states, which means that the errors in the band gap and defect levels within the gap do not cancel in total-energy differences. These error are not easily corrected (see [this paper on correction methods](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.78.235104){:target="_blank"}).  

The best approach to handling this error is to use HSE total energies for each of the systems. This HSE total energy should only be taken from one k-point because the defect level will have some artificial dispersion due to the finite supercell size. We use the $\Gamma$-point because, in general, the error at gamma is considered to be the smallest. If that is found to not be the case, another k-point should be used to further minimize the error.

If you are using a $\Gamma$-only k-point grid, the eigenvalues from the HSE level should be used as they are most accurate. For larger grids, however, the PBE eigenvalues can be used because the differences between bands within the conduction and valence bands are relatively accurate and HSE calculations at all k-points may not be feasible. 

## VASP

### Recommended order

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

{% include note.html content="The [`Export`][code_export] can currently only export the wave functions from the standard version of VASP. The energies can be extracted from the standard and $\Gamma$-only versions. Assumptions made in the code are not compatible with the noncollinear version." %}

For HSE calculations, start from a preconverged PBE calculation for speed. It is fastest to do HSE calculations with the $\Gamma$-only version of VASP, but keep in mind the current limitations of the [`Export`][code_export] code. If you do use the $\Gamma$-only version to get accurate HSE energies, make sure that your preconverged PBE calculation was run using the $\Gamma$-only version as well, otherwise you will get an error saying `plane wave coefficients changed`.

### Parameters

#### [`EDIFF`](https://www.vasp.at/wiki/index.php/EDIFF){:target="_blank"}

For the relaxations, it is okay to use a looser `EDIFF` like `1e-5`, but we use the wave functions from the SCF calculations, so the convergence criteria must be tighter for those. We have used `1e-8` in the past. To ensure that your convergence is tight enough, run the [`TME`][code_tme] code with the same system input for both `PC` and `SD` and a small set of bands used for the initial state and final state. Check the `|<f|i>|^2` column in the `allElecOverlap.*.*` file to ensure that the bands are orthonormal. For the same initial and final band, the overlap should be close to 1. For different initial and final bands, the overlap should be close to zero (`~1e-10` is usually good).

#### [`ISYM`](https://www.vasp.at/wiki/index.php/ISYM){:target="_blank"}

Our code is not currently implemented with symmetrization, so symmetry must be completely turned off with `ISYM=-1`.

#### [`ISPIN`](https://www.vasp.at/wiki/index.php/ISPIN){:target="_blank"}

Make sure that the proper choice of `ISPIN` is used for each of your calculations. The [`Export`][code_export] and [`TME`][code_tme] codes handle spin polarization automatically, so it is not required to use the same setting for every system. 

#### [`NBANDS`](https://www.vasp.at/wiki/index.php/NBANDS){:target="_blank"}

From the VASP wiki:
> In the electronic minimization calculations, empty states do not contribute to the total energy, however, empty states are required to achieve a better convergence. 

The default setting in VASP already includes extra bands above the filled states to get better convergence; however, we need accurate calculations of both occupied and unoccupied states when considering electron capture. For electron capture, ensure that `NBANDS` is set so that there are plenty empty bands _above the empty states that you are using for your initial electron states._ The default value from VASP is `NELECT/2+NIONS/2`. You should at least increase above the default by the number of unoccupied states that you need.

#### [`NCORE`](https://www.vasp.at/wiki/index.php/NCORE){:target="_blank"}

`NCORE` sets the number of bands per group. It is related to `NPAR`, which is the number of bands that are treated in parallel. The recommended settings on the wiki are `NPAR = sqrt(number-of-cores)` or `NCORE = cores-per-node`. The relationship between `NPAR` and `NCORE` is `NCORE = number-of-cores/KPAR/NPAR`. 

We have found that, for the system sizes we typically work with (~100-atom supercells), `NCORE=~16` is good for PBE and `NCORE=~4` is good for HSE. Ideally, `NCORE` should be a factor of `cores-per-node`, since this reduces communication between nodes. A lower (higher) `NCORE` is slower (faster) but utilizes less (more) memory.

`NBANDS` and `NCORE` should line up as follows:
  * Take the number of cores and divide by `NCORE` to get the number of band groups
  * The number of bands must be evenly divisible by the number of band groups
  * If these values don't line up, VASP will automatically adjust the numbers of bands, which can lead to issues if you are setting the occupations manually using `FERWE/FERDO`
  
#### [`KPAR`](https://www.vasp.at/wiki/index.php/KPAR){:target="_blank"}

`KPAR` should be an integer divisor of the total number of cores. The number of processors working on a group of k-points is `number-of-cores/KPAR`. Try not to let the number of processors working on a group split across nodes (i.e., make `number-of-cores/KPAR` an integer divisor of the number of cores per node). Also, try to have the same number of k-points in each group (i.e., make `KPAR` an integer divisor of the number of unique k-points) or at least balance the load as best as possible. As an example, with 10 k-points you could have 5 groups working on 2 k-points each, but you wouldn't want 8 groups because then 6 groups would be idle while the two groups with 2 k-points work on their second point. 

In VASP, the k-points are split up first, then bands. To be safe, you can set the bands and processors based on `NCORE`, then set `KPAR` based on number of k-points and multiply the number of nodes by `KPAR`

{% include links.html %}

## Export

Once the VASP calculations are complete, the information needs to be exported to be fed into the next calculations in line. For [zeroth-order capture][zeroth-order-capture], you only need the wave functions from the perfect crystal and the final charge state in the relaxed initial positions. 

{% include note.html content="The [`Export`][code_export] assumes that the VASP calculation was SCF and not relax, otherwise total energy may be extracted incorrectly; However, you should always do an SCF calculation after relaxation anyways if you need anything other than the final positions." %}

To export the wave functions, you should use an input file like
```fortran
&inputParams
  VASPDir = 'path-to-VASP-output'                       ! default './'
  exportDir = 'path-to-put-exported-files'              ! default './Export'
  gammaOnly = .false.
  energiesOnly = .false.
/
```

For systems that you only need to extract the energies from (total energy and eigenvalues), use `energiesOnly = .true.`. The extracted energies get fed into the [`EnergyTabulator`][calculations_energy-tabulator] to tabulate the various energies needed for different parts of the theory.
