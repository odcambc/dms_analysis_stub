# DMS Analysis Stub

## Description
This repository contains some files and structure that can be used to initialize a new DMS analysis project. It is intended to be used as a template for new projects.

This workflow contains the starting framework used in the Fraser and Coyote-Maestes labs for analyzing DMS data. We will attempt to keep it up to date with the current versions of our basic analyses.

## Table of Contents
- [DMS Analysis Stub](#dms-analysis-stub)
  - [Description](#description)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
  - [Usage](#usage)
    - [File Structure](#file-structure)
  - [Citations](#citations)
  - [License](#license)
  - [Contributing](#contributing)

## Installation
Clone the repository into a new directory:
```
git clone https://github.com/odcambc/dms_analysis_stub.git
```
Alternatively, fork the repository and clone your fork. This will allow you to make changes to the stub and use your own version. Please consider contributing your changes back to the main repository, as well!

After cloning, you will need to initialize the repository. If you are using rstudio, you can do this by opening the project file `dms_analysis_stub.Rproj`, then initializing the environment with the comment `renv::restore()`.

Otherwise, you can initialize the repository by running the following command from the root directory of the repository:

```
Rscript -e 'renv::restore()'
```

To start analyzing and plotting, simply your data files to the `data` folder and edit the analysis files in the `analysis` directory. Use the file structure below as a guide for where to put your files.

## Usage
### File Structure
The file structure can be modified as needed, but the default structure is as follows:
```
.
├── analysis
│   ├── functions
│   │   └── dms_analysis_utilities.r
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

## Citations
This workflow, in some version, has been used in the following manuscripts:

- [Rao et al, 2023](https://www.biorxiv.org/content/10.1101/2023.10.24.562292v1)
- [Estevam et al, 2023](https://www.biorxiv.org/content/10.1101/2023.08.03.551866v1)
- [Yee et al, 2023](https://www.biorxiv.org/content/10.1101/2023.06.06.543963v1)
- [Macdonald et al, 2022](https://www.biorxiv.org/content/10.1101/2022.07.26.501589v1)

If you found this useful, please let us know!

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing
Pull requests and issues are welcome and appreciated.