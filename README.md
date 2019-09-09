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

### Examples

This basic example shows how to extract the genomic features (_i.e._, the genes, with their upstream and downstream regions) from a TxDb object

```
library("Mimir")

library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")

query = extractGenomicFeatures(TxDb = TxDb.Dmelanogaster.UCSC.dm3.ensGene)
```

> Extracted 10562 genes
> 
> Extracted 10559 upstream regions
> 
> Extracted 10562 downstream regions

Now that the regions have been extracted, the signal across these features can be calcualted

```
library("GenomicAlignments")
library("pasillaBamSubset")

# Locate and read the alignemnt file
fl1 <- untreated1_chr4()
subject = readGAlignments(fl1)

profile = profileGenomicFeatures(genomicRegions=query, sampleObject=subject, TxDb=TxDb.Dmelanogaster.UCSC.dm3.ensGene)
```

> Upstream regions on + strand: 5261
>
> Upstream regions on - strand: 5298
> 
> Gene_body regions on + strand: 5264
> 
> Gene_body regions on - strand: 5298
> 
> Downstream regions on + strand: 5264
> 
> Downstream regions on - strand: 5298

Finally, the profile obtained can be visualised

```
library("ggplot2")

ggplot(profile, aes(x=bin, y=Mean, colour=region_id)) + 
   geom_line() +
   geom_vline(xintercept=c(50.5, 150.5), linetype="dashed", colour="grey30") +
   scale_x_continuous("Relative position",
        breaks=c(1, 50.5, 150.5, 250), label=c("-500", "TSS", "TES", "1000")) +
   scale_y_continuous("Average normalised signal") +
   coord_cartesian(xlim=c(0, 250)) +
   theme_bw() +
   theme(legend.position=c(0.9, 0.8), legend.background=element_blank()) +
   guides(colour=guide_legend(title=""))
```

![](man/figures/example_profile.png)

## Authors

* **Federico Agostini**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
