# DMS Analysis Stub

## Description
This repository contains files and a structure that can be used to initialize a new DMS analysis project. It is intended to be used as a template for new projects.

  This workflow contains the starting framework used in the [Fraser](https://fraserlab.com/) and [Coyote-Maestas](https://www.wcoyotelab.com/) labs at UCSF for analyzing DMS data. We will attempt to keep it up to date with the current versions of our basic analyses.

## Table of Contents
- [DMS Analysis Stub](#dms-analysis-stub)
  - [Description](#description)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
  - [Usage](#usage)
    - [File Structure](#file-structure)
    - [Custom functions](#custom-functions)
  - [Citations](#citations)
  - [License](#license)
  - [Contributing](#contributing)

## Installation
Clone the repository into a new directory:
```
git clone https://github.com/odcambc/dms_analysis_stub.git
```
Alternatively, fork the repository and clone your fork. This will allow you to make changes to the stub and use your own version. Please consider contributing your changes back to the main repository, as well!

After cloning, you will need to initialize the repository. If you are using rstudio, you can do this by opening the project file `dms_analysis_stub.Rproj`, then initializing the environment with the command `renv::restore()`.

Otherwise, you can initialize the repository by running the following command from the root directory of the repository:

```
Rscript -e 'renv::restore()'
```

To start analyzing and plotting, simply add your data files to the `data` folder and edit the analysis files in the `analysis` directory. Use the file structure below as a guide for where to put your files.

## Usage
### File Structure
The file structure can be modified as needed, but the default structure is as follows:
```
.
├── analysis
│   ├── functions
│   │   ├── dms_analysis_utilities.r 
│   │   ├── map_pdb.r
│   │   ├── mkh_heatmap_function.R
│   │   ├── plot_heatmap.r
│   │   └── plot_scatter.r
│   ├── dms_cutoffs.Rmd
│   └── dms_plotting.Rmd
├── data
│   ├── scores
│   │   └── (Put the starting score files in here)
│   ├── sequences
│   │   └── (Put the WT sequence fasta in here)
│   └── structures
│       └── (Put any PDB files in here)
├── dms_analysis_stub.Rproj
├── LICENSE
├── README.md
├── renv
│   └── ...
└── results
    ├── plots
    │   ├── heatmaps
    │   │   └── (Output heatmaps will be put in here)
    │   └── scatterplots
    │       └── (Output scatterplots will be put in here)
    ├── scores
    │   └── (Output score files will be put in here)
    └── structures
        └── (Output structures will be put in here)
```

### Custom functions

The `analysis/functions` directory contains some custom functions that are used in the analysis. Here is a brief description of some of the utilities:

`dms_analysis_utilities.r` contains general functions for parsing DMS data and setting up aesthetics. Color palettes and
  some useful variant orderings for plotting are also defined here.

`map_pdb.r` contains a function for mapping positionally summarized DMS data onto a PDB structure. This function uses the `bio3d` package.

`plot_heatmap.r` contains a function for plotting a heatmap of DMS data. This function is fairly customizable, but the defaults are hopefully reasonable. Some useful variant orderings for plotting are defined in `dms_analysis_utilities.R`.

`mkh_heatmap_function.R` contains an alternative heatmap plotting function. Note that the syntax is _not_ identical to `plot_heatmap.r`!

`plot_scatter.r` contains a function for plotting scatterplots and marginal densities of DMS data. This function also tries
to determine some rough significance thresholds for scores using
the synonymous variant effect distribution. See the comments for more details.

## Citations
This workflow, in some version, has been used in the following manuscripts:

- [Rao et al, 2023](https://www.biorxiv.org/content/10.1101/2023.10.24.562292v1)
- [Estevam et al, 2023](https://www.biorxiv.org/content/10.1101/2023.08.03.551866v1)
- [Yee et al, 2023](https://www.biorxiv.org/content/10.1101/2023.06.06.543963v1)
- [Macdonald et al, 2022](https://www.biorxiv.org/content/10.1101/2022.07.26.501589v1)

If you found this useful, please let us know!

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.

## Contributing
Pull requests and issues are welcome and appreciated.
