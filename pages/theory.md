---
title: "Theory"
keywords: theory
tags: [theory]
sidebar: home_sidebar
permalink: theory.html
summary: This page gives a summary of the theory developed in the Barmparis and GaN papers.
---

## Background

There are two kinds of capture: equilibrium and nonequilibrium. In equilibrium capture, the carrier is near the defect and the rates of capture and emission are equal. This case is covered by the capture theory dating back to the 1950s with [Huang and Rhys](https://royalsocietypublishing.org/doi/abs/10.1098/rspa.1950.0184){:target="_blank"}, where the carrier is perturbed by electron-phonon coupling. 

However, most experimental conditions will have some kind of current, making the system nonequilibrium. For the nonequilibrium case, we must consider a carrier coming far away from the defect and being perturbed by the presence of the defect. Perturbation theory cannot be directly used for carrier for capture because the defect state cannot be constructed from a finite basis set of perfect-crystal Bloch states. 

The theory behind the capture formalism (originally developed in the [Barmparis paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.214111){:target="_blank"}) is based on defining a basis set where the initial band state is "decoupled" from the final defect band state. The wave functions are defect states with the coupling matrix element between the initial and final states removed. The coupling matrix element is used as a perturbation, and the full Hamiltonian is the full defect system.The transition rate is derived using Fermi's golden rule. The wave functions start as many-body but the form of Fermi's golden rule allows the orthonormality of the single-particle orbitals to be exploited to simplify (assuming the form of a Slater determinant). 

The wave functions discussed above are only the electronic states. The electronic and vibrational states are decoupled using the <a href="#" data-toggle="tooltip" data-original-title="{{site.data.glossary.boa}}">Born-Oppenheimer approximation.</a> The phonons come from the defect system. The phonon overlaps are transformed using a Fourier-like transform to a time-domain integral. Feynman's path integral is then used to simplify the integrand. 

Starting from the "uncoupled" basis set and transforming to usable wave functions (defect and pristine- crystal states), we end up with two terms. One term is first-order in the atomic displacements. This is exactly the term from electron-phonon coupling from the previous theory. Additionally, however, we end up with a zeroth-order (in the atomic displacements) term not considered by any other paper. This term represents the nonequilibrium perturbation on a carrier coming far from the defect. 

The first-order term was first calculated using density functional theory by [Alkauskas et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.075202){:target="_blank"} for nonradiative capture in 2014. However, that approach and all previous approaches to the phonon integral have used single-mode or effective-mode approximations because the number of possible phonon configurations blows up as more modes are allowed to participate. We have developed a time-domain integration form that allows explicit consideration of all of the phonon modes in the crystal. 

Below, I discuss the form of the zeroth- and first-order terms and how they relate to our code.  

## Zeroth-order

The final form of the transition rate for the zeroth-order is 

$$\Gamma_i^{(0)} = \frac{2}{\hbar^2}g \left|M_{\text{e}}^{\text{BO}} ( \{ \mathbf{R}_k^{(i)} \} ) \right|^2 \text{ Re} \int_0^{\infty} G^{(0)}(t) e^{i E_{if}^{\text{DFT}} t/\hbar - \gamma t} \,dt$$

where $g$ is a degeneracy factor for the summing over (nearly) degenerate final electronic states. I will explain the main components of this equation and how it relates to the code in the sections below.

{% include note.html content="The derivation of the zeroth-order term assumes that the carrier is coming from far away from the defect and gets perturbed by the coupling with the defect as it approaches. This treatment is only valid for nonequilibrium capture. In principle, most experiments will be nonequilibrium because they will have some current; however, when comparing to truly equilibrium conditions, the zeroth-order term should not be present. The first-order term should be present in either case because the perturbation from electron-phonon coupling always exists." %}

### Electronic matrix element and energy differences

$M_{\text{e}}^{\text{BO}}( \\{ \mathbf{R}_k^{(i)} \\} )$ is the zeroth-order matrix element:

$$M_{\text{e}}^{\text{BO}} ( \{ \mathbf{R}_k^{(i)} \} ) = \langle\phi_f| \phi_i^0 \rangle \left[ E_{f}(\{ {\bf R}_k^{(i)}\}) - E_{i}(\{ {\bf R}_k^{(i)}\}) \right]$$

This matrix element is calculated in the [`TME`][code_tme] code. $\|\phi_f \rangle$ is the single-particle final defect state taken from the [`Export`][code_export] of the SCF calculation on the __final charge state in the initial positions__. $\|\phi_i^0 \rangle$ is the single-particle initial perfect-crystal Bloch state taken from the [`Export`][code_export] of the SCF calculation on the pristine supercell. 

The energy difference in $M_{\text{e}}^{\text{BO}}$ is the many-body electronic-only energy difference from just the carrier changing state. We use the total energy difference between the two charge states in the initial positions:

$$W_{f}( \{ {\bf R}_k^{(i)}\}) - W_{i}( \{ {\bf R}_k^{(i)} \} ) \approx W_{f}( \{ {\bf R}_k^{(i)} \} ) - \left[ W_{\text{ref}}( \{ {\bf R}_k^{(i)} \} ) + \varepsilon_{i} - \varepsilon_{\text{ref}} \right].$$

### Time-domain integral and kernel

The time-domain integral comes from a Fourier-like transform of the delta function. The kernel of the integration, $G^{(0)}(t)$, is derived with help from Feynman's path integral and a transformation to $q$ space. It is defined as

$$G^{(0)}(t) = \exp \left[ \sum_j S_j(\bar{n}_j + 1) e^{i\omega_j t}+S_j \bar{n}_j e^{-i\omega_j t}-S_j(2\bar{n}_j+1)\right],$$

where

$$S_j = \frac{\omega_j}{2\hbar}(\Delta q_j)^2$$

is the Huang-Rhys factor and

$$\bar{n}_j = \frac{1}{e^{\beta \hbar \omega_j}-1}.$$

The Huang-Rhys factor includes $\Delta q_j$, which is the projection of the displacement vector between the relaxed positions of the two charge states onto the phonon eigenvectors of mode $j$. The $\Delta q_j$, $\omega_j$, and $S_j$ are tabulated in the [`PhononPP`][calculations_phononpp] code. The final transition rate is then calculated by the [`LSF`][code_lsf] code, which primarily does the time-domain integral.


{% include links.html %}
