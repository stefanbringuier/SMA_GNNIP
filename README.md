<p align="center">
<br>
<br>
<a href="https://github.com/stefanbringuier/SMA_Phonons_GNNIP/actions/workflows/build.yml">
<img src="https://github.com/stefanbringuier/SMA_Phonons_GNNIP/actions/workflows/build.yml/badge.svg?branch=main" alt="Article status"/>
</a>
<a href="https://github.com/stefanbringuier/SMA_Phonons_GNNIP/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/article-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/stefanbringuier/SMA_Phonons_GNNIP/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
</p>

# Assessing Phonon Properties of Shape Memory Alloys using Graph Neural Network Potentials

**Authors:** Stefan Bringuier

**Affiliation:** Independent Researcher, San Diego, CA

## Abstract
Shape Memory Alloys (SMAs) such as NiTi are pivotal in vibrational dampening,robotic arms, prothestic hands, and mechanical vales as a result of their unique displacive phase transformation properties. However,assessing their phonon characteristics, which provide insight into the displacive phase transformation, by using conventional Density Functional Theory (DFT) methods is computationally challenging. In this paper we assess the use of recent Graph Neural Network (GNN) potentials as a viable alternative. These GNN potentials have leveraged large DFT-computed databases for training and provide a good approximation to thepotential energy surfaces which enable inference on energies and forces in lattice and atomic dynamic calculations. We compare GNN potentials with the Embedded Atom Method (EAM) in predicting phonon properties of NiTi. We find that in general the EAM potentials are slightly more performant than the GN in that they appear to more accurately describe equation of state, phonon dispersions, and elastic constants. However, the MACE potential shows significant promise as an exploritory potential for designing and characterizing more complex SMAs. Our findings reveal GNN potentials as a promising avenue for advancing SMA research and fine-tuning of pre-trained GNN models with more relevant DFT data (i.e., force-constants) may permit better outcomes.


## Reproduce
If you want to reproduce the raw data, plots, and manuscript from scratch:

1. Install `showyourwork`
2. Clone this repo.
3. From the command-line run `showyourwork clean` and then `showyourwork build`

This may take quite some time (3+ hours) depending on your resources. Note that all the GNN potentials run by default on the cpus, and this code base has no flags for switching it to `cuda`. However, this should be relatively simple by just modifying the [Calculators.py](src/scripts/Calculators.py) script. The phonon calculations use considerable memory and you'll need at least 32 GB of RAM.

### Dataset
The data for this manuscript is stored in the file [NiTi_Structures.json](src/data/NiTi_Structures.json). This is a human readable [ASE style JSON](https://wiki.fysik.dtu.dk/ase/ase/db/db.html) file that can be either used with a standard JSON library or ASE. Within the database are ASE `Atoms` type structures and then additional datafields with regard to the GNN model used, spacgroup of the structure, equation of state parameters, phonon bandstructure/dos, and elastic constants.

### Plotting
Any file in [src/scripts](src/scripts)  with a variation of the word "Plot" generates some plot. Not all plots are used in the manuscript, for specifics see the [Snakefile](Snakefile) or [showyourwork.yml](showyourwork.yml) files.

### Manuscript
A single `LaTeX` files is used, [ms.tex](src/tex/ms.tex) and is fairly boilerplate other than the written text. The rendered figures are saved in [src/tex/figures](src/tex/figures) and are in `.png` format. Some formatting with regard to placement of subfigures is done in the `ms.tex` document within the `\begin{figure} ... \end{figure}` environment. So the figures in the manuscript pdf may be combinations of figures in [src/tex/figures](src/tex/figures)

### Github Actions
Because multiple conda environments are used and the packages are considerable in size, the disk space used during the Github actions is very large. Therefore the hack provided by <https://github.com/marketplace/actions/maximize-build-disk-space> needs to be used.

### Acknowledgement
An open source scientific article created using the [showyourwork](https://github.com/showyourwork/showyourwork) workflow.
