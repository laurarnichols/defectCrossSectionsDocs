---
title: "Tutorial: $\\text{Si}_{\\text{V}}\\text{H}_3$"
keywords: tutorial silicon example
tags: []
sidebar: home_sidebar
permalink: tutorial.html
summary: This page includes a tutorial on running an electron capture calculation from start to finish using a triply-hydrogenated silicon vacancy as an example.
---

## Problem definition

This tutorial will consider the positive-to-neutral transition of a triply-hydrogenated silicon vacancy. We use the [WZP](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.105501){:target="_blank"} method to accomplish our desired charge states without affecting the neutrality of the supercell. The ground state of this defect is neutral with a singly-occupied defect. To get the positive charge state, we excite the electron out of the defect level into the conduction band. 

## VASP and Export

0. Make sure the executables are set up on your machine (see [Getting Started](getting-started))
1. Do a volume relaxation on the silicon unit cell. The conventional unit cell that we start with is
```
Si
1.000000
   5.430700000  0.000000000  0.000000000
   0.000000000  5.430700000  0.000000000
   0.000000000  0.000000000  5.430700000
   Si
    8
Direct
0.000000000 0.000000000 0.000000000
0.000000000 0.500000000 0.500000000
0.500000000 0.000000000 0.500000000
0.500000000 0.500000000 0.000000000
0.250000000 0.250000000 0.250000000
0.250000000 0.750000000 0.750000000
0.750000000 0.250000000 0.750000000
0.750000000 0.750000000 0.250000000
```
You can also start from the primitive cell as in [this VASP tutorial](https://www.vasp.at/wiki/index.php/Cd_Si_volume_relaxation){:target="_blank"}. You should end up with a lattice parameter of $5.4684 \,\overset{\lower.5em\circ}{\mathrm{A}}$.
2. Scale up to your supercell (we used 4x4x4, i.e. 512 atoms). 
3. Relax the inner DOF of the perfect crystal. For this tutorial, just use the gamma-point, but make sure to use the standard version of VASP for any calculations you plan to export the wave functions from (pristine and final charge state at initial positions). The `INCAR` used for the pristine relaxation was 
```fortran
! How to initialize charge density and wave functions
ISTART = 0        ! 0 to start a new job, 1 for continuing and reading in wf
ICHARG = 2        ! 2 to start a new job, 1 for continuing and reading in chg
INIWAV = 1        ! chooses initially random wf, numerically best
!
! Electronic control flags
ENCUT  = 400      ! cutoff energy
EDIFF  = 1E-6     ! convergence threshold
PREC   = Normal   ! numerical precision
ALGO   = Normal   ! Davidson algorithm, usually best choice
NELMIN = 5        ! min electronic steps
NELM   = 50       ! max electronic steps
!
ISPIN  = 1        ! Spin polarized calculation: 1 = No, 2 = Yes
!
! Parallelization flags
NCORE = 11        ! Band parallelization, 16 is good for PBE, 4 is good for HSE
KPAR  = 1         ! Set to number of k-point blocks (remember to multiply nodes/procs in jobscript)
!
! Ionic control flags
EDIFFG = -1E-2    ! Convergence for ionic relaxation
NSW = 300           ! Number of ionic steps per relaxation job
IBRION = 2       ! Relaxation algorithm, -1 is off, 2 is conjugate gradient
ISIF = 2          ! What to relax. 2 does just the internal DOF, 3 does the lattice vectors too
!
! Electronic broadening
ISMEAR = 0        ! 0 is Gaussian broadening
SIGMA  = 0.01     ! amount of broadening
!
LWAVE  = .TRUE.   ! write out wavefunctions
LCHARG = .TRUE.   ! write out chargedensity
LVTOT  = .FALSE.  ! write out full potential
LVHAR  = .FALSE.  ! write out electrostatic potential
```
_The `!` characters on the blank lines are needed for the above code block to be formatted correctly, but they serve no purpose in the actual file._ __Only the gamma point should be used for all of the VASP calculations in this example.__
4. Follow with an SCF calculation with a tighter `EDIFF=1e-8` and more bands. We used around 1500, but make sure that the number of bands that you choose fits the number of processes you are running on (see the [`NBANDS` and `NCORE` sections on the VASP and Export page](vasp-and-export.html#nbands)). 
5. In the center of the supercell, introduce the defect. In this case, remove a silicon atom close to the center and add 3 hydrogen atoms along the previously-existing bonds. One silicon atom will retain a dangling bond, which is the defect level. 
6. Relax the inner DOF in the ground state, then run the required SCF calculation to fix the charge density for the last relaxation step.
7. Use the converged `WAVECAR` and `CHGCAR` as the starting point for the postive-charge-state relaxation. To get the positive charge state, use
```fortran
ISMEAR = -2
FERWE = 1022*1 1 0 0.5 535*0
FERDO = 1022*1 1 0 0.5 535*0
```
`FERWE` sets the spin-up occupations and `FERDO` sets the spin-down occupations. The numbers given are for a choice of 1560 bands. Make sure that your number of bands aligns with your [parallelization criteria](vasp-and-export.html#ncore), otherwise VASP may automatically adjust it, messing up the manual occupations. I like to have the band extrema and the defect levels written explicitly for easier visualization in different files, but that is not required. Note that the excited electron is split across the spin channels to remove any artificial spin-orbit coupling with the defect.
8. Use the relaxed positive-charge state positions and run a ground-state SCF calculation with the same parameters as the pristine SCF calculation.
9. [`Export`](vasp-and-export.html#export) the pristine supercell calculation and the ground-state-defect/positive-state-positions to get the wave functions for the zeroth-order term (you only need to set `VASPDir` and `exportDir` for these calculations).

To get the total energies, also do HSE SCF calculations for the relaxed positive charge state, the relaxed ground state, and the ground state at the positive relaxed positions. The gamma-only version of VASP can be used since we only need the energies. [`Export`](calculations_export) the HSE results with the option `energiesOnly = .true.` and `gammaOnly =.true.`, if applicable. 

## EnergyTabulator

Once all of the results are exported, run [`EnergyTabulator`](calculations_energy-tabulator) to calculate all of the required energies. Make sure to point `exportDirEigs` to the ground-state-HSE `Export` output. Since we are using HSE total energies and eigenvalues, we don't need any energy corrections. The reference band should be set to where the carrier is in our calculations. In this case, the reference carrier is in the lowest conduction band at the $\Gamma$-point:
```fortran
refBand = 1025
```
We are considering electron capture, so we should also set
```fortran
captured = .true.
elecCarrier = .true.
```
The initial and final band ranges should be set based on the initial and final states that you want to consider. For capture, there should be only one final state (the defect level), so we should use
```fortran
iBandFinit = 1024
iBandFfinal = 1024
```
Include as many initial bands as you want to plot. We used 
```fortran
iBandIinit = 1025
iBandIfinal = 1088
```
which is about a 1 eV range to plot. This calculation is really fast ($<$ 1 min.), so you can run on a single node. 

## Phonons and PhononPP

Use some software to get the ground-state phonon eigenvectors and frequencies. We use [Phonopy](https://phonopy.github.io/phonopy/){:target="_blank"} (see [our documentation](phonons) on how to run a simple phonon calculation), but any software can be used as long as the result is the same YAML format. 

Once you have the phonons, run the [`PhononPP`](calculations_phononpp) code to get all of the shifted positions needed for the first-order matrix elements and calculate the Huang-Rhys factor among other things needed as input to the [`LSF`][code_lsf] code. _You will want to run this in a clean directory because a shifted `POSCAR` will be output for each phonon mode._ We used `shift = 0.01` $\overset{\lower.5em\circ}{\mathrm{A}}$.

## Zeroth-order

Getting the [zeroth-order][zeroth-order-capture] transition rate is pretty straightforward. First, run the [`TME`][code_tme] code with `order=0`, with `exportDirSD` pointing to the [`Export`][code_export] output for the defect in the ground charge state and initial positions and `exportDirPC` pointing to the `Export` output for the perfect crystal. The same band ranges used in the `EnergyTabulator` code should be used here, and you should point to the energy tables. 

Once that `TME` is complete, run [`LSF`][code_lsf] with `order=0` and the same band ranges. Make sure to also point to the energy tables with `energyTableDir`, the `TME` output with `matrixElementDir`, and the `Sj.out` file generated by the `PhononPP` code using `SjInput`. The key parameters in the `LSF` input file are the temperature, time step, and smearing parameters (see the [`LSF`][code_lsf] documentation for more information on these parameters). Use
```fortran
 temperature          = 300
 dt                   = 2.1d-9
 hbarGamma            = 1.0    ! meV
 smearingExpTolerance = 1d-4
```
It is also important to pick out the proper spin channel. For this problem, choose the first spin channel with `iSpin=1`.

## First-order

The [first-order term][first-order-capture] is more involved and much longer (total around 4 days). 

### TME

To calculate the first-order matrix element, the following calculations must be done for each shifted file for the different phonon modes:
* VASP SCF calculation
* [`Export`][code_export]
* [`TME`][code_tme]

A separate directory must be created for each mode and the shifted `POSCAR` must be copied into the directory along with the required VASP files. __Make sure that each folder has padding zeros to make the number 4-digits wide, otherwise there will be issues with reading the files in the `LSF` code.__ To speed things up, it also helps to copy the `WAVECAR` and `CHGCAR` from the unshifted calculation. It is best to set up a bash script to do this rather than manually copying all of the files. If you include copying the `WAVECAR` and `CHGCAR`, this will take a long time, so it is best to submit the script to the queue to allow it to run. 

#### Shifted overlaps

For increased accuracy, it is recommended to first run [`TME`][code_tme] using the nondisplaced system as both inputs to get numbers that represent numerical "orthonormality." This serves as the zero for the displaced overlap; subtracting the baseline will improve the numerical accuracy of the matrix elements. Calculate these overlaps first, then use `subtractBaseline = .true.` and include the path to the basline overlaps in `baselineDir`.

The results for the paper used a version of the [`TME`][code_TME] code that did not have the ability to subtract out a baseline for the overlaps (now implemented), so I had to manually subtract out the baseline using my helper script. I left the script I originally used to shift the matrix elements in case it might be helpful as a starting point to do some post-processing in the future because the first-order matrix elements are very costly to repeat. I would highly suggest using [ChatGPT](https://chat.openai.com/){:target="_blank"} to help generate any scripts you need if you are not familiar with bash scripting, but beware that you will still need to check the results from ChatGPT as it is not always completely accurate. 


#### Helper script outline

I utilized a [bash script](https://github.com/laurarnichols/defectCrossSections/blob/daa793dac9d4ea30617e40ddfb8eb04ad5cdfce8/Tools/firstOrderHelper.scr){:target="_blank"} (used on NERSC's Perlmutter) with different functions for different processes (more details below):
* `submitVASP`
    * Set up the displacement directories
    * Copy needed files for VASP calculation into folder
    * Generate a VASP run script with the job labeled by the mode number
    * Submit VASP
* `clean`
    * Change into a displacement directory and delete all files and the `export` folder
* `submitAll`
    * Do everything as in `submitVASP`
    * Also copy the input file for `Export` and generate the input file for `TME`
    * Create and submit a script to run VASP, `Export`, and `TME`, with the latter two only running if the previous calculation submitted without errors
* `submitExportAndTME`
    * Copy and generate `Export` and `TME` input files
    * Create and submit a script to run `Export` then `TME` if `Export` completes properly
* `cleanForResubmit`
    * Remove all files and folders from a displacement directory
* `getShiftedOverlaps`
    * Copy the header from the original `TME` output file
    * Get $\delta q_j$ displacement from the `dq.txt` file
    * Paste together into a temporary file the data needed to calculate a corrected matrix element:
        * Overlap squared from displaced calculation
        * Overlap squared from non-displaced calculation (baseline)
        * First-order eigenvalue difference
        * Displacement for this mode (repeated for every band)
    * Recalculate the matrix element and put into a temporary file
    * Get all of the data from the original `TME` output file but replace the matrix element column with the corrected matrix element from the temporary file

The loop at the bottom of the file allows the above functions to be run for multiple modes (make sure to use leading zeros to match the file names). When copying over wave functions and charge densities to submit VASP calculations for multiple modes, make sure to submit the script to the queue.

{% include note.html content="Make sure to be considerate of the machine you are working on. Do not abuse the debug queue or submit too many calculations at once (I maxed out around 250). Do not take up too much storage space. Make sure to delete the VASP and `Export` files when you are sure `TME` is complete. These calculations are huge, and you do not want to get throttled or lose your access because of bad practices." %}

#### Ideal workflow

The ideal workflow would be to use `submitAll` to submit all calculations and just track them as they finish. You can tell the `Export` and `TME` codes are finished by looking for a line like
```
************ VASP Export complete! (     35.01 secs) ************
```
for `Export` or
```
-----------------------------------------------------------------
 Total time needed:                              63.02 secs.
```
for `TME`. However, there are a few things that require active monitoring.

#### Confirm calculations were successful

Sometimes something will be wrong with the machine that will cause calculations to hang or die, so you need to double-check that each of the calculations finished correctly. I did this using a few different commands to search through the output files. 
* `General timing` in the VASP `OUTCAR`
* `Export complete` in the `Export` standard out file
* `Total time` in the `TME` standard out file

The presence of those lines means that each calculation completed successfully. It is not feasible to look through every file individually if you are doing many modes, so it is best to use the command line to process the files. 

As an example, in the ideal workflow, VASP, `Export`, and `TME` would all be run one after the other. I would start with checking 100 files at a time to see if `TME` completed.
```bash
grep "Total time" disp-00*/TME.out | wc -l
```
The `grep` command searches for the phrase in a given path. Specifying `disp-00*` means that the first 99 directories (only 99 for `00` because there is no mode `0`) will be searched. The `wc -l` gets the line count from the output, representing in how many files the line was found. If the number represents all of the files searche by the wildcard, you can mark all of those modes as complete. If not, search by tens (e.g., `disp-001*`) to narrow down which files did not complete properly. If a group of ten is not all complete, you can remove the `| wc -l` part to list all of the results and show which modes completed properly. 

__I would highly recommend using a spreadsheet or something similar to track the status of the matrix elements for all of the modes.__

### LSF

Once all of the first-order matrix elements are complete, run [`LSF`][code_lsf] with `order=1` and the same parameters as in the zeroth-order term. `matrixElementDir` should be set to the folder used for the matrix-element output in each subfolder (usually `TME`). `MjBaseDir` should be set to the main folder that all of the displaced-calculation folders are in. `prefix` should be set to the prefix of each subdirectory (usually `disp-`).

## Group velocity

The results from the LSF code give you the transition rate, but the transition rate can be converted to a cross section using 

$$\sigma_i = \frac{\Gamma_i \Omega}{v_g^{(i)}},$$

where $\Gamma_i$ is the transition rate, $\Omega$ is the supercell volume, and $v_g^{(i)}$ is the group velocity for the perfect-crystal band state $i$. 

To accurately get the group velocity, it is best to utilize the equation

$$v_{gx}^{(i)} = \frac{1}{\hbar} \frac{\partial E_i}{\partial k_x},$$

for each direction $x,y,z$, then use

$$||v_g|| =  \sqrt{v_{gx}^2 + v_{gy}^2 + v_{gz}^2}$$

to combine the components. To practically calculate this, we must get the perfect-crystal energy bands at k-points slightly offset from the base k-point and calculate the slope of the bands. You only need two points to determine a slope, but there are often band crossings. Often, three points (center and $\pm$) is enough to connect straight lines, but sometimes the bands are curved and using straight lines would give incorrect results. To be safe, 4 points should be used in each coordinate direction so that it is clear how the points should be connected and how the slope should be calculated. 

As an example, consider the `KPOINTS` file centered at the `0.000000 0.000000 0.000000` k-point (in direct coordinates):
```
Explicit k-points
10
Direct
 0.000000  0.000000  0.000000  1
-0.010000  0.000000  0.000000  0
 0.010000  0.000000  0.000000  0
 0.020000  0.000000  0.000000  0
 0.000000 -0.010000  0.000000  0
 0.000000  0.010000  0.000000  0
 0.000000  0.020000  0.000000  0
 0.000000  0.000000 -0.010000  0
 0.000000  0.000000  0.010000  0
 0.000000  0.000000  0.020000  0

```
The main k-point is the only one with any weight to speed up the calculations. Each direction contains one negatively-displaced k-point, two positively-displaced k-points, and the center k-point. Run an SCF calculation on the perfect crystal supercell with the given set of k-points. 

Once that calculation is complete, you need to gather the eigenvalues for plotting so that you can manually determine how the bands should line up. You can do this manually by looking at the `EIGENVAL` file, or you can use the [`Export`][code_export] code. To do this, use the following options:
```fortran
energiesOnly = .true.
groupForGroupVelocity = .true.
nDispkPerCoord = 3
pattern = '-0.01 0.01 0.02'
```
The `Export` code assumes that the base k-point is first in the list and that the displacement pattern for each coordinate is consistent. It does not matter what order you do the displacement in, as long as the pattern is consistent for each coordinate and it matches the `pattern` variable in the input file. `pattern` should have `nDispkPerCoord` numbers in a string. The velocities will be sorted and output from negative to positive (including the base k-point) for each coordinate in files like `groupedEigenvaluesX.isp.ik`. 

{% include note.html content="There are symmetry considerations for calculating the group velocity at the $\Gamma$ point. I do not know how this was done, but Guanzhi is supposed to be writing up documentation on this process." %}

I created a [jupyter notebook](https://github.com/laurarnichols/CrossSectionCalculations/blob/4dc5181bff428e8244c92048b9ff0726454181d0/SiVH3/GaNPaper/VASPandExport/pristine/groupVelocity/GroupVelocity.ipynb){:target="_blank"} to help speed up the process of lining up bands and calculating the slope. You will need to manually swap bands at the different k-points so that each band is a smooth curve. You will need to do this for each direction $x,y,z$. You may also need to swap some bands across the coordinates so that the total velocity $v_g^{(i)}$ is the same for all bands within a degeneracy group.

## Smearing and plotting

Finally, the results must be smeared based on the group velocity to produce a smooth curve as a function of energy. I created a [jupyter notebook](https://github.com/laurarnichols/CrossSectionCalculations/blob/4dc5181bff428e8244c92048b9ff0726454181d0/SiVH3/GaNPaper/results/SmearCrossSectionWithVg.ipynb){:target="_blank"} to plot the smeared cross section and also explore the different elements that go into the results.

{% include links.html %}