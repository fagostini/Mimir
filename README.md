# Mimir package

R package to create metadata profiles for CLIP and other types of NGS data.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisities

A pre-installed version of [R](https://www.r-project.org/) (>= 3.6.0; it has not been tested on any previous version, but it might work), with a working version of [Bioconductor](https://www.bioconductor.org/).

Additional dependencies will be installed via BiocManager:install() 

### Install

__[Option 1]__ Use _install_github("author/package")_.

Install the _devtools_ package and load it

```
install.packages("devtools")
library(devtools)
```

Install the [fagostini/Mimir](https://github.com/fagostini/Mimir) package with

```
install_github("fagostini/Mimir", type = "source", repos = BiocManager::repositories(), dependencies = TRUE)
```

__[Option 2]__ Use _install.packages()_

Download the source release from [here](https://github.com/fagostini/Mimir/releases) and install it with

```
install.packages(path_to_file, type = "source", repos = BiocManager::repositories(), dependencies = TRUE)
```

## Authors

* **Federico Agostini**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
