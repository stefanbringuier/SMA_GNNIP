<p align="center">
<br>
<br>
<!---
<a href="https://github.com/stefanbringuier/SMA_Phonons_GNNIP/actions/workflows/build.yml">
<img src="https://github.com/stefanbringuier/SMA_Phonons_GNNIP/actions/workflows/build.yml/badge.svg?branch=main" alt="Disabled"/>
</a>
-->
<a href="https://github.com/stefanbringuier/SMA_Phonons_GNNIP/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/preprint_status-work%20in%20progress-orange.svg?style=flat" alt="WIP"/>
</a>
<a href="https://github.com/stefanbringuier/SMA_Phonons_GNNIP/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/build-disabled-white.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/stefanbringuier/SMA_Phonons_GNNIP/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/preprint-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/stefanbringuier/SMA_Phonons_GNNIP/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/preprint-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
</p>

<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css">
<div class="alert alert-warning" role="alert" style="font-family: Arial, sans-serif; font-size: 0.8em;">
    <strong>⚠️ This research is <em>in progress</em>, therefore the results/analysis are still being completed. Please also excuse any spelling and grammar mistakes within the preprint. ⚠️</strong>
</div>



## Suitability of Graph Neural Network for Shape Memory Alloys

**Authors:** Stefan Bringuier

**Affiliation:** Independent Researcher, San Diego, CA

**Contact:** [stefanbringuier@gmail.com](mailto:stefanbringuier@gmail.com)

### Abstract

> Shape Memory Alloys (SMAs), such as NiTi, are crucial in vibration dampening, robotic arms, prosthetic hands, and mechanical valves due to their unique displacive phase transformation properties. However, assessing their phonon characteristics using conventional Density Functional Theory (DFT) methods is computationally expensive. This paper evaluates the use of recent Graph Neural Network (GNN) potentials as a viable alternative. These GNN potentials, trained on extensive DFT-computed databases, approximate potential energy surfaces, enabling inference on energies and forces in lattice and atomic dynamics calculations. We compare GNN potentials with the Embedded Atom Method (EAM) and Modified EAM (MEAM) potentials in predicting phonon properties of NiTi and PtTi. Generally, EAM/MEAM potentials appear to slightly outperform GNNs phonon, providing more accurate descriptions of the equation of state, phonon dispersions, and elastic constants. However, the MACE potential shows significant promise as an exploratory tool for designing and characterizing more complex SMAs. Our findings suggest GNN potentials as a promising initial avenue for advancing SMA research, particularly for ternary and quaternary systems, with improved outcomes possible by fine-tuning pre-trained GNN models with more relevant DFT data (e.g., force-constants).


## Reproduce Manuscript

If you want to reproduce the raw data, plots, and manuscript from scratch:

1. Install `showyourwork`
2. Clone this repo.
3. From the command-line run `showyourwork clean` and then `showyourwork build`

This may take quite some time (3+ hours) depending on your resources. The conda environments take some time to configure. If you can use the `--cores NUMBER_CPUS --conda-frontend=mamba` options.  Note that all the GNN potentials run by default on the `cpus`, and this code base has no flags for switching it to `cuda`. However, this should be relatively simple by just modifying the [Calculators.py](src/scripts/Calculators.py) script. The phonon calculations use considerable memory and you'll need at least 32 GB of RAM.

### Dataset

The data for this manuscript is stored in the file [SMA_Results.json](src/data/SMA_Results.json)[^1]. This is a human readable [ASE style JSON](https://wiki.fysik.dtu.dk/ase/ase/db/db.html) file that can be either used with a standard JSON library or ASE. Within the database are ASE `Atoms` type structures and then additional datafields with regard to the GNN model used, spacgroup of the structure, equation of state parameters, phonon bandstructure/dos, and elastic constants.

### Plotting

Any file in [src/scripts](src/scripts)  with a variation of the word "Plot" generates some plot. Not all plots are used in the manuscript, for specifics see the [Snakefile](Snakefile) or [showyourwork.yml](showyourwork.yml) files.

### Manuscript

A single `LaTeX` files is used, [ms.tex](src/tex/ms.tex) and is fairly boilerplate other than the written text. The rendered figures are saved in [src/tex/figures](src/tex/figures) and are in `.png` format. Some formatting with regard to placement of subfigures is done in the `ms.tex` document within the `\begin{figure} ... \end{figure}` environment. So the figures in the manuscript pdf may be combinations of figures in [src/tex/figures](src/tex/figures)

### Github Actions
> Disabled for now due to disk size limitations for Github actions.
>
Because multiple conda environments are used and the packages are considerable in size, the disk space used during the Github actions is very large. Therefore the hack provided by <https://github.com/marketplace/actions/maximize-build-disk-space> needs to be used.

### Collaborations

If your a researcher who wants to collaborate please message me.

#### License
The code in this repository is open sourced under the [MIT license](LICENSE).

#### Acknowledgement

This reproducible workflow and manuscript is enabled by [showyourwork](https://github.com/showyourwork/showyourwork). I would also like to thank [Dr. Venkateswara Manga](https://mse.engineering.arizona.edu/faculty-staff/faculty/venkateswara-manga) for fruitful discussions.

#### Footnotes
[^1]: This is will be generated when you run `showyourwork build`.
