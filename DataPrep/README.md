
"cmwf_csv.sh" is the processing pipeline that takes as input HTS
*.fastq files as well as various reference files, and produces as
output various summary files and *.vcf files. Please see block
comments within the script for more information.

This pipeline has several dependencies, which are, in effect,
described by the Dockerfile.

"VCF_mutation_picker.0.6.py" is invoked as the last step of the
pipeline. It produces a per-gene summary of results.


Creating an image file
======================

If you have docker, singularity and suitable permissions[^1], it is straightforward to
create an image with all pipeline dependencies installed.

In the following we
use "cmwf" as a tag/name for the docker and singularity images. Check
to make sure this doesn't collide with an image already known to your
docker setup.

1. Build the initial docker image\
`docker build -t cmwf - < Dockerfile`

2. Instantiate a container (this will place in a shell---just exit from it)\
`docker run -it cmwf`

3. Find the size of the image using\
`docker image`

4. Find the container id using\
`docker ps -a`

5. Create an image file (replace 1500 if need be with a number in MB that is slightly more than the size reported in step (3))\
`singularity image.create -s 1500 /your/image/directory/cmwf.img`

6. Transform from the image from docker to singularity (note the "sudo")
`docker export <container id from step (4)> | sudo $(which singularity) image import /your/image/directory/cmwf.img`

Once the final image is created, feel free to remove the intermediate docker container and image ("docker rm" and "docker rmi" commands).


Running the pipeline via an image file
======================================

In these examples we assume you have a directory setup something like:

`$FREYA_ROOT/FREYA
$FREYA_ROOT/data
$FREYA_ROOT/image
$FREYA_ROOT/results`

Where the first is a clone of this repo, the second holds the phenotype csv file and a subdirectory with the fastq files, the third your image file, and the fourth will hold the results.

There are two common execution modes:

1. Using the resources of a single node. In this example the pipeline script will run entirely within a singularity session. The environment variable `DISBATCH_SSH_NODELIST` is set to indicate that the subtasks should be run on the local node.

```bash
DISBATCH_SSH_NODELIST=localhost:$(nproc) singularity exec -B ${FREYA_ROOT} ${FREYA_ROOT}/image/cmwf.img bash ${FREYA_ROOT}/FREYA/cmwf_csv.sh ${FREYA_ROOT}/data/phenotype.csv ${FREYA_ROOT}/data/fastqs ${FREYA_ROOT}/results```

This would be a good starting point with a small test data set, but will take a while with a large one.

2. Using the resources of multiple nodes, here via SLURM. In this case the driver script is run via a SLURM submission script, but each of the pipeline tasks is executed within a singularity session (specified by the environment variable `DB_TASK_PREFIX`) is set to indicate that the subtasks should be run on the local node. The submission script invokes the pipeline driver script using:

```bash
DB_TASK_PREFIX="singularity exec -B ${FREYA_ROOT} ${FREYA_ROOT}/image/cmwf.img " singularity exec -B ${FREYA_ROOT} ${FREYA_ROOT}/image/cmwf.img bash ${FREYA_ROOT}/FREYA/cmwf_csv.sh ${FREYA_ROOT}/data/phenotype.csv ${FREYA_ROOT}/data/fastqs ${FREYA_ROOT}/results```

Because the driver script will be run in your normal (non-singularity) environment via a batch submission script, in this case you must have [disBatch](https://github.com/flatironinstitute/disBatch) installed and in your PATH.



[^1]: If push comes to shove, you could use could use a system like VirtualBox to create a virtual machine with docker, singularity and the appropriate permissions. This VM is needed only to build the image, not to execute it.

