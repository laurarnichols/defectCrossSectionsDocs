---
title: "Tutorial: $\\text{Si}_{\\text{V}}\\text{H}_3$"
keywords: tutorial silicon example
tags: []
sidebar: home_sidebar
permalink: tutorial.html
summary: This page includes a tutorial on running an electron capture calculation from start to finish using a triply-hydrogenated silicon vacancy as an example.
---

## Problem definition

This tutorial will consider the positive-to-neutral transition. We use the [WZP](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.105501){:target="_blank"} method to accomplish our desired charge states without affecting the neutrality of the supercell. The ground state of this defect is neutral with a singly-occupied defect. To get the positive charge state, we excite the electron out of the defect level into the conduction band. 

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
4. Follow with an SCF calculation with a tighter `EDIFF=1e-8` and more bands. We used around 1500, but make sure that the number of bands that you choose fits the number of processes you are running on (see the [`NBANDS` and `NCORE` sections on the VASP page](vasp.html#nbands)). 
5. In the center of the supercell, introduce the defect. In this case, remove a silicon atom close to the center and add 3 hydrogen atoms along the previously-existing bonds. One silicon atom will retain a dangling bond, which is the defect level. 
6. Relax the inner DOF in the ground state, then run the required SCF calculation to fix the charge density for the last relaxation step.
7. Use the converged `WAVECAR` and `CHGCAR` as the starting point for the postive-charge-state relaxation. To get the positive charge state, use
```fortran
ISMEAR = -2
FERWE = 1022*1 1 0 1 535*0
FERDO = 1022*1 1 0 0 535*0
```
`FERWE` sets the spin-up occupations and `FERDO` sets the spin-down occupations. The numbers given are for a choice of 1560 bands. Make sure that your number of bands aligns with your [parallelization criteria](vasp.html#ncore), otherwise VASP may automatically adjust it, messing up the manual occupations. I like to have the band extrema and the defect levels written explicitly for easier visualization in different files, but that is not required.
8. Use the relaxed positive-charge state positions and run a ground-state SCF calculation with the same parameters as the pristine SCF calculation.
9. [`Export`](calculations_export) the pristine supercell calculation and the ground-state-defect/positive-state-positions (you only need to set `VASPDir` and `exportDir` for these calculations).

To get the total energies, also do HSE SCF calculations for the relaxed positive charge state, the relaxed ground state, and the ground state at the positive relaxed positions. The gamma-only version of VASP can be used since we only need the energies. [`Export`](calculations_export) the HSE results with the option `energiesOnly = .true.` and `gammaOnly =.true.`, if applicable. 

## EnergyTabulator

Once all of the results are exported, run [`EnergyTabulator`](calculations_energy-tabulator) to calculate all of the required energies. Make sure to point `exportDirEigs` to the ground-state-HSE `Export` output. Since we are using HSE total energies and eigenvalues, we don't need any energy corrections. The reference band should be set to where the carrier is in our calculations. In this case, the reference carrier is in the CBM:
```fortran
refBand = 1025
```
We are considering electron capture, so we should also set
```fortran
CBMorVBMBand = 1025
```
The initial and final band ranges should be set based on the initial and final states that you want to consider. For capture, there should be only one final state (the defect level), so we should use
```fortran
iBandFinit = 1024
iBandIinit = 1024
```
Include as many initial bands as you want to plot. We used 
```fortran
iBandIinit = 1025
iBandIfinal = 1088
```
which is about a 1 eV range to plot. This calculation is really fast ($<$ 1 min.), so you can run on a single node. 

## Phonons and PhononPP

Use some software to get the ground-state phonon eigenvectors and frequencies. We use [Phonopy](https://phonopy.github.io/phonopy/){:target="_blank"} (see [our documentation](phonons) on how to run a simple phonon calculation), but any software can be used as long as the result is the same YAML format. 

Once you have the phonons, run the [`PhononPP`](calculations_phononpp) code to get all of the shifted positions needed for the first-order matrix elements and calculate the Huang-Rhys factor among other things needed as input to the [`LSF`][calculations_lsf] code. _You will want to run this in a clean directory because a shifted `POSCAR` will be output for each phonon mode._ We used `shift = 0.01` $\overset{\lower.5em\circ}{\mathrm{A}}$.

## Zeroth-order

Getting the zeroth-order transition rate is pretty straightforward. First, run the [`TME`][calculations_tme] code with `order=0`, with `exportDirSD` pointing to the [`Export`][calculations_export] output for the defect in the ground charge state and initial positions and `exportDirPC` pointing to the `Export` output for the perfect crystal. The same band ranges used in the `EnergyTabulator` code should be used here, and you should point to the energy tables. 

Once that `TME` is complete, run [`LSF`][calculations_lsf] with `order=0` and the same band ranges. Make sure to also point to the energy tables with `energyTableDir`, the `TME` output with `matrixElementDir`, and the `Sj.out` file generated by the `PhononPP` code using `SjInput`. The key parameters in the `LSF` input file are the temperature, time step, and smearing parameters (see the [`LSF`][calculations_lsf] documentation for more information on these parameters). Use
```fortran
 temperature          = 300
 dt                   = 2.1d-9
 hbarGamma            = 1.0    ! meV
 smearingExpTolerance = 1d-4
```
It is also important to pick out the proper spin channel. For this problem, choose the spin-down channel with `iSpin=1`.

## First-order

The first-order term is more involved. To calculate the first-order matrix element, the following calculations must be done for each shifted file for the different phonon modes:
* VASP SCF calculation
* [`Export`][calculations_export]
* [`TME`][calculations_tme]

A separate directory must be created for each mode and the shifted `POSCAR` must be copied into the directory along with the required VASP files. __Make sure that each folder has padding zeros to make the number 4-digits wide, otherwise there will be issues with reading the files in the `LSF` code.__ To speed things up, it also helps to copy the `WAVECAR` and `CHGCAR` from the unshifted calculation. It is best to set up a bash script to do this rather than manually copying all of the files. If you include copying the `WAVECAR` and `CHGCAR`, this will take a long time, so it is best to submit the script to the queue to allow it to run. Here is an example PBS script of how that can be done:
```bash
#!/bin/bash
# Move POSCAR-xxx files from `PhononPP` into individual folders and
# create a run script for each. You must have INCAR, KPOINTS, and
# POTCAR files in current directory. If you want to include the 
# WAVECAR and CHGCAR as well, make sure those are present in the 
# current directory or you include the proper path.
#
# To run (assuming you call this script `createShiftedDirs.scr`):
# chmod +x createShiftedDirs.scr
# qsub createShiftedDirs.scr
#
## Required PBS Directives --------------------------------------
#PBS -A <grant-id>
#PBS -q standard
#PBS -l select=3:ncpus=44:mpiprocs=44
#PBS -l walltime=1:00:00
#PBS -l place=scatter:excl

## Optional Directives ------------------------------------
#PBS -N SiVH3_submit
#PBS -j oe
#PBS -M <email>
#PBS -m be

# Define parameters
BASE_DIR="${WORKDIR}/path/to/your/folder"

N=132


cd ${BASE_DIR}


createDirCopy() {
  local i=$1
  # Make folder for each displacement
  if [ ! -d "disp-$i" ]; then
    mkdir disp-$i

    # Move the displaced POSCAR files to individual
    # folders as POSCAR
    mv ph_POSCAR_$i disp-$i/POSCAR

    # Copy in other input files for SCF calculation, Export, and TME
    cp INCAR KPOINTS POTCAR WAVECAR CHGCAR SiVH3_export.in TME.in disp-$i

    # Generate a script to run all 3 in this folder
    cat > disp-$i/run_VASP-Export-TME.scr << EOF
#!/bin/bash
## Required PBS Directives --------------------------------------
#PBS -A <grant-id>
#PBS -q standard
#PBS -l select=6:ncpus=44:mpiprocs=44
#PBS -l walltime=1:00:00
#PBS -l place=scatter:excl

## Optional Directives ------------------------------------
#PBS -N ${i}SiVH3_FOME
#PBS -j oe

# Define parameters
BASE_DIR="\${WORKDIR}/path/to/your/folder/disp-${i}"

## Execution Block -----------------------------------------------
# Environment Setup
# cd to base directory in /work
cd \${BASE_DIR}

# Launching  ----------------------------------------------------
module unload cray-libsci

echo "Running VASP." &&
date +"%T" &&
aprun -n 264 ~/vasp.5.4.4/bin/vasp_std &&
echo "VASP exited clean. Running Export." &&
date +"%T" &&
aprun -n 264 ~/defectCrossSections/bin/Export_VASP.x -nk 1 -nb 4 < SiVH3_export.in > SiVH3_export.out &&
echo "Export exited clean. Submitting TME." &&
date +"%T" &&
aprun -n 264 ~/defectCrossSections/bin/TME.x -nk 1 < TME.in > TME.out &&
date +"%T"
EOF
  fi
}

for i in {0001..1539}
do
  createDirCopy $i &

  if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
    wait -n
  fi
done
wait
```

This script will utilize multiple processes to copy the directories faster. You can edit this script or create others for submitting the scripts, checking their status, copying the output, etc. This is a huge calculation (total of 4ish days), and it requires checking on the calculations regularly. I would highly suggest using [ChatGPT](https://chat.openai.com/){:target="_blank"} to help generate any scripts if you are not familiar with bash scripting, but beware that you will still need to check the results from ChatGPT as it is not always completely accurate. 

Some tips for running these calculations:
* Track which ones are complete in a spreadsheet or something similar.
* You can run these calculations on multiple machines if you have access to speed things up.
* If you run on multiple machines, make sure to track which calculations are on which machine (they will all need to be gathered on one machine at some point).
* Use `wc -l` on the output files to check that the calculations finished (use a calculation you know completed as reference).

Once all of the first-order matrix elements are complete, run [`LSF`][calculations_lsf] with `order=1` and the same parameters as in the zeroth-order term. `matrixElementDir` should be set to the folder used for the matrix-element output in each subfolder (usually `TME`). `MjBaseDir` should be set to the main folder that all of the displaced-calculation folders are in. `prefix` should be set to the prefix of each subdirectory (usually `disp-`).

## Group velocity

The results from the LSF code give you the transition rate, but the transition rate can be converted to a cross section using 

$$\sigma = \frac{\Gamma \Omega}{v_g},$$

where $\Gamma$ is the transition rate, $\Omega$ is the supercell volume, and $v_g$ is the group velocity. 

To accurately get the group velocity, it is best to utilize the equation

$$v_{gi} = \frac{1}{\hbar} \frac{\partial E}{\partial k_i},$$

where $i$ is the $x,y,z$ index, then use

$$||v_g|| =  \sqrt{v_{gx}^2 + v_{gy}^2 + v_{gz}^2}$$

to combine the components. To practically calculate this, we must get the perfect-crystal energy bands at k-points slightly offset from the base k-point and calculate the slope of the bands. You only need two points to determine a slope, but there are often band crossings. Often, three points (center and $\pm$) is enough to connect straight lines, but sometimes the bands are curved and using straight lines would give incorrect results. To be safe, 5 points should be used in each coordinate direction so that it is clear how the points should be connected and how the slope should be calculated. 

As an example, consider the `KPOINTS` file centered at the `0.000000 0.000000 0.000000` k-point (in direct coordinates):
```
Explicit k-points
13
Direct
 0.000000  0.000000  0.000000  1
-0.020000  0.000000  0.000000  0
-0.010000  0.000000  0.000000  0
 0.010000  0.000000  0.000000  0
 0.020000  0.000000  0.000000  0
 0.000000 -0.020000  0.000000  0
 0.000000 -0.010000  0.000000  0
 0.000000  0.010000  0.000000  0
 0.000000  0.020000  0.000000  0
 0.000000  0.000000 -0.020000  0
 0.000000  0.000000 -0.010000  0
 0.000000  0.000000  0.010000  0
 0.000000  0.000000  0.020000  0

```
The main k-point is the only one with any weight to speed up the calculations. Each direction contains two negatively-displaced k-points, two positively-displaced k-points, and the center k-point. Run an SCF calculation on the perfect crystal supercell with the given set of k-points. 

Once that calculation is complete, you need to gather the eigenvalues for plotting so that you can manually determine how the bands should line up. You can do this manually by looking at the `EIGENVAL` file, or you can use the [`Export`][calculations_export] code. To do this, use the following options:
```fortran
energiesOnly = .true.
groupForGroupVelocity = .true.
nDispkPerCoord = 4
pattern = '-0.02 -0.01 0.01 0.02'
```
The `Export` code assumes that the base k-point is first in the list and that the displacement pattern for each coordinate is consistent. It does not matter what order you do the displacement in, as long as the pattern is consistent for each coordinate and it matches the `pattern` variable in the input file. `pattern` should have `nDispkPerCoord` numbers in a string. The velocities will be sorted and output from negative to positive (including the base k-point) for each coordinate in files like `groupedEigenvaluesX.isp.ik`. 

{% include note.html content="There are symmetry considerations for calculating the group velocity at the gamma point. I do not know how this was done, but Guanzhi is supposed to be writing up documentation on this process." %}

{% include links.html %}