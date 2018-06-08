
[![Binder](http://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/flatironinstitute/FREYA/master)

# FREYA


This is a public repository for canine breast cancer analysis.

There are two main sub-pipelines in this analysis. The first, in DataPrep, will process raw data to create expression and mutation calls. The second will take that data (or you can provide your own processed data) and perform a series of analyses, including comparing it to human breast cancer. 

The DataAnalysis part of the repository is designed to use 
[`repo2docker`](https://repo2docker.readthedocs.io/en/latest/)
to build a `docker` container which includes
the code for running the analysis and all needed dependencies.
The `repo2docker` tool requires the `docker` infrastructure
and Python 3 to run.  See the 
[installation instructions](https://repo2docker.readthedocs.io/en/latest/install.html).

Directions for running and installing both DataPrep and DataAnalysis are located in their respective folders. 


## Todo (temp):

Fix notebook errors.

Output artifact (like diagrams) should go into a mounted persistent location.

Make repo public after vetting.

Provide mechanisms for uploading and downloading data that will work in Binder environment.
