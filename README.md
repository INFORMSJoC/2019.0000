[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Data and R scripts for the Case Study of Unmanned Aerial Vehicle Information Collection Missions with Uncertain Characteristics

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The R scripts and data in this repository are a snapshot of the scripts and data that were used in the case study reported on in the paper
[Unmanned Aerial Vehicle Information Collection Missions with Uncertain Characteristics](https://doi.org/) by M.D. Moskal, E. Dasdemir and R. Batta.


## Cite
To cite this material, please cite the [paper](https://doi.org/) using its DOI and the software itself, using the following DOI.

[![DOI](https://zenodo.org/badge/285853815.svg)](https://zenodo.org/badge/latestdoi/285853815)

Below is the BibTex for citing this version of the code.

```
@article{CacheTest,
  author =        {Michael D. Moskal, Erdi Dasdemir, Rajan Batta},
  publisher =     {INFORMS Journal on Computing},
  year =          {2022},
  doi =           {},
  url =           {https://github.com/INFORMSJoC/},
}  
```

## Description

The goal of this repository is to share data and R scripts of the case study of our paper.  Our motivation is to present our code and results in a reproducible way and facilitate the coding effort of those who want to run further experiments or improve our model.

## Repository Structure
We have a quite simple repository structure.

- [data and scripts](data_and_scripts) folder contains two subfolders for the civilian and military case application. Each subfolder contains the input data in text files and the corresponding R script. There is also an additional spreadsheet (.xlsx file) in the main folder with raw data to make it easier to understand the data. Note that Each R script reads the data from the input files in the same location. Please set the working directory accordingly at the first lines of the scripts.

## Data
We adapt the data set provided by [Dasdemir et al. (2022)](https://doi.org/10.1016/j.cor.2022.105882) to create the problem instances for civilian and military applications. This data set includes the elevation and coordinate data of 58 fourteeners, and the water capacity, elevation, and coordinate data of 41 water reservoirs from the State of Colorado. Please see [data and scripts](data_and_scripts) and [spreadsheet file](data_and_scripts/colorado-data-for-uav-case-study.xlsx) for the data.

## R scripts
R scripts employ [Gurobi](www.gurobi.com) solver to solve the mathematical optimization model. You should first install Gurobi on your computer and then the Gurobi R package to your R insallation.

- Please visit [Gurobi web page](https://www.gurobi.com/) for installing Gurobi and obtaining a license.
-  Please visit [Gurobi R Installation Guide](https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html) and [Installing the R Package](https://www.gurobi.com/documentation/9.5/quickstart_windows/r_ins_the_r_package.html) pages to use Gurobi on R.

## Support
For support in using the scripts and data you can reach Erdi Dasdemir (erdidasdemir@gmail.com).
