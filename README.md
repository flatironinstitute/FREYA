
[![Binder](http://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/flatironinstitute/FREYA/master)

# FREYA


This is a public repository for FREYA, a computational pipeline for cross-species cancer analyses. 

The FREYA analytic framework calculates genetic diversity within a cancer cohort, extracts progression-related patterns of expression, calls functional variants, identifies intrinsic human molecular tumor subtypes, and compares them to human breast cancer biology. It helps streamline future canine mammary tumor (or other tissue/species) analyses and provides a comprehensive suite of tools that encompass conventional human analyses and new dog-centric approaches, seamlessly integrating the two.

There are two main sub-pipelines in this analysis. The first, DataPrep, processes raw sequencing data to create expression and mutation calls. The second, DataAnalysis, takes that data (or you can provide your own processed data) and performs a series of analyses, including comparisons to human breast cancer.

Directions for running and installing [DataPrep](https://github.com/flatironinstitute/FREYA/blob/master/DataPrep/README.md) and [DataAnalysis](https://github.com/flatironinstitute/FREYA/blob/master/DataAnalysis/README.md) are located in their respective folders. There is also an overview of the framework and its requirements at [freya.flatironinstitute.org](http://freya.flatironinstitute.org/). 


