---
title: "Theory"
keywords: theory
tags: [theory]
sidebar: home_sidebar
permalink: theory.html
summary: This page gives a summary of the theory developed in the Barmparis and GaN papers.
---

The theory behind the capture formalism (originally developed in the [Barmparis paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.214111){:target="_blank"}) is based on defining a basis set where the initial band state is decoupled from the final defect band state. This definition is based on the physical problem of a carrier coming far from the defect and getting captured. The problem is treated through perturbation theory. The unperturbed set of wave functions is the decoupled basis set, and the full wave functions are the defect wave functions. 

The transition rate is derived using Fermi's golden rule. The wave functions start as many-body but the form of Fermi's golden rule allows the orthonormality of the single-particle orbitals to be exploited to simplify (assuming the form of a Slater determinant). These wave functions are only the electronic states. The electronic and vibrational states are decoupled using the <a href="#" data-toggle="tooltip" data-original-title="{{site.data.glossary.boa}}">Born-Oppenheimer approximation.</a> The phonons come from the defect system. The phonon overlaps are transformed using a Fourier-like transform to a time-domain integral. Feynman's path integral is then used to simplify the integrand. 

## Zeroth-order

The final form of the transition rate for the zeroth-order is 

$$\Gamma_i^{(0)} = \frac{2}{\hbar^2} \left|M_{\text{e}}^{\text{BO}} ( \{ \mathbf{R}_k^{(i)} \} ) \right|^2 \text{ Re} \int_0^{\infty} G^{(0)}(t) e^{i E_{if}^{\text{DFT}} t/\hbar - \gamma t} \,dt.$$

I will explain each component of this equation and how it relates to the code in the sections below.

{% include note.html content="The derivation of the zeroth-order term assumes that the carrier is coming from far away from the defect and gets perturbed by the coupling with the defect as it approaches. This treatment is only valid for nonequilibrium capture. In principle, most experiments will be nonequilibrium because they will have some current; however, when comparing to truly equilibrium conditions, the zeroth-order term should not be present. The first-order term should be present in either case because the perturbation from electron-phonon coupling always exists." %}

### Electronic matrix element and energy differences

$M_{\text{e}}^{\text{BO}}( \\{ \mathbf{R}_k^{(i)} \\} )$ is the zeroth-order matrix element:

$$M_{\text{e}}^{\text{BO}} ( \{ \mathbf{R}_k^{(i)} \} ) = \langle\phi_f| \psi_i^0 \rangle \left[ E_{f}^{\text{DFT}}(\{R_k^{(i)}\}) - E_{i}^{\text{DFT}}(\{R_k^{(i)}\}) \right]$$

This matrix element is calculated in the [`TME`][calculations_tme] code. $\|\phi_f \rangle$ is the single-particle final defect state taken from the [`Export`][calculations_export] of the SCF calculation on the __final charge state in the initial positions__. $\|\psi_i^0 \rangle$ is the single-particle initial perfect-crystal Bloch state taken from the [`Export`][calculations_export] of the SCF calculation on the pristine supercell. 

The energy difference in $M_{\text{e}}^{\text{BO}}$ is the many-body electronic-only energy difference from just the carrier changing state. It is based on the total DFT energy difference between the defect system in the different charge states but in the initial relaxed positions ($\\{R_k^{(i)}\\}$). The other energy ($E_{if}^{\text{DFT}}$) in the exponential is defined as 

$$E_{if}^{\text{DFT}} \equiv E_{f}^{\text{DFT}}(\{R_k^{(f)}\}) - E_{i}^{\text{DFT}}(\{R_k^{(i)}\}),$$

which is the total energy difference between the two relaxed charge states. Both energies are calculated in the [`EnergyTabulator`][calculations_energy-tabulator] code. More details on the exact form of the energies are given on that documentation page.

### Time-domain integral and kernel

The time-domain integral comes from a Fourier-like transform of the delta function. The kernel of the integration, $G^{(0)}(t)$, is derived with help from Feynman's path integral and a transformation to $q$ space. It is defined as

$$G^{(0)}(t) = \exp \left[ \sum_j S_j(\bar{n}_j + 1) e^{i\omega_j t}+S_j \bar{n}_j e^{-i\omega_j t}-S_j(2\bar{n}_j+1)\right],$$

where

$$S_j = \frac{\omega_j}{2\hbar}(\Delta q_j)^2$$

is the Huang-Rhys factor and

$$\bar{n}_j = \frac{1}{e^{\beta \hbar \omega_j}-1}.$$

The Huang-Rhys factor includes $\Delta q_j$, which is the projection of the displacement vector between the relaxed positions of the two charge states onto the phonon eigenvectors of mode $j$. The $\Delta q_j$, $\omega_j$, and $S_j$ are tabulated in the [`PhononPP`][calculations_phononpp] code. The final transition rate is then calculated by the [`LSF`][calculations_lsf] code, which primarily does the time-domain integral.


{% include links.html %}
