==========
PELE basis
==========

Our Protein Energy Landscape Exploration (PELE) tool combines protein structure prediction algorithms and Metropolis Monte Carlo techniques
to efficiently tackle tasks like predicting binding sites, docking pose refinement or modelling exit path of a ligand.

PELE algorithm
--------------

Each simulation consists of several steps executing the following algorithm:

**1. Perturbation.** Localised perturbation of the ligand (if present), involving random translation and rotation,
followed by simple side chain relocation to avoid clashes. Additionally, the protein is minimized by driving alpha
carbons to new positions resulting from a small displacement in a low frequency anisotropic normal mode (ANM).

**2. Relaxation.** Optimization of side chains in proximity to the ligand as well as those whose energy changed the
most during ANM, using a rotamer library with a resolution of 10°. This is followed by a global minimization with
Truncated Newton minimizer.

**3. Acceptance.** The new structure is accepted (defining a new minimum) or rejected based on the Metropolis criterion.

.. image:: https://nostrumbiodiscovery.github.io/pele_docs/_images/pele_scheme.png
  :width: 400
  :align: center

Further reading
---------------

Over the years, numerous publications have been written about the methodology and applications of PELE itself, as well
as further improvements, such as AdaptivePELE or FragPELE.

* `PELE: Protein Energy Landscape Exploration. A Novel Monte Carlo Based Technique <https://pubs.acs.org/doi/abs/10.1021/ct0501811>`_ by Kenneth W. Borrelli, Andreas Vitalis, Raul Alcantara, and Victor Guallar

* `Adaptive simulations, towards interactive protein-ligand modeling <https://www.nature.com/articles/s41598-017-08445-5>`_ by Daniel Lecina, Joan F. Gilabert, and Victor Guallar

* `FragPELE: Dynamic Ligand Growing within a Binding Site. A Novel Tool for Hit-To-Lead Drug Design <https://pubs.acs.org/doi/10.1021/acs.jcim.9b00938>`_ by Carles Perez, Daniel Soler, Robert Soliva, and Victor Guallar
