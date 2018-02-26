
# FREYA


This is a public repository for canine breast cancer analysis.

The repository is designed to use 
[`repo2docker`](https://repo2docker.readthedocs.io/en/latest/)
to build a `docker` container which includes
the code for running the analysis and all needed dependencies.
The `repo2docker` tool requires the `docker` infrastructure
and Python 3 to run.  See the 
[installation instructions](https://repo2docker.readthedocs.io/en/latest/install.html).

[Temporary note: eventually the repo should have the option to use Binder
with no local install required.]

To build the docker image from the command line

```bash
jupyter-repo2docker --no-run \
    --image-name canine_bc \
    https://github.com/flatironinstitute/FREYA
```

The analysis requires source data which is not included
in the repository to protect privacy and intellectual property.
The required data must be mounted to the container or provided
in some other way. We provide synthetic sample data in TODO.

To run the image with a mounted data directory located at
`/Path/to/your/data`
invoke `docker run` with the directory mounted like so:

```bash
docker run -it --rm -p 8888:8888 \
   -v /Path/to/your/data:/home/rstudio/notebooks/user_data:z \
   canine_bc:latest
```

Changes made in the data directory by the image will persist in the mounted directory.

When the docker instance starts up the log output looks like this:

```
[I 19:34:40.443 NotebookApp] Writing notebook server cookie secret to /home/rstudio/.local/share/jupyter/runtime/notebook_cookie_secret
[I 19:34:40.935 NotebookApp] Serving notebooks from local directory: /home/rstudio
[I 19:34:40.935 NotebookApp] 0 active kernels
[I 19:34:40.935 NotebookApp] The Jupyter Notebook is running at:
[I 19:34:40.935 NotebookApp] http://0.0.0.0:8888/?token=7baca1829f0665df26adf9943c6f04279647171eb3ffd12c
[I 19:34:40.935 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[W 19:34:40.935 NotebookApp] No web browser found: could not locate runnable browser.
[C 19:34:40.936 NotebookApp] 
    
    Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://0.0.0.0:8888/?token=7baca1829f0665df26adf9943c6f04279647171eb3ffd12c
```

Paste the `http:...` URL into your browser including your token (which will be different
from the one shown above) to connect to the Jupyter server.  

## Using the server in the browser

The Jupyter server provides a browser based interface to the directories
in the container and the notebooks and other files in those directories.
Please see the
[notebook basics notebook](http://nbviewer.jupyter.org/github/jupyter/notebook/blob/master/docs/source/examples/Notebook/Notebook%20Basics.ipynb)
for an introduction to using the Jupyter server web interface.

## Shutting down the server

To terminate the Jupyter server and the Docker instance run type `Control-C` twice in the
console window where you started the instance.


## Connection trouble?

If you have trouble connecting to the jupyter server it may be because some other process
is responding to port 8888 on your system.  To find out whether processes are connected on
that port use

```bash
lsof -ti:8888 
```

If you are sure the processes are not needed you
can use `kill -9` to terminate any such processes.

```bash
$ lsof -ti:8888 
1510
$ kill -9 1510
```


## For developers: mounting the notebooks for persistence

To build the docker image from a local directory including provisional changes
that have not been checked in to the github repository use 
`repo2docker` with the path to the directory
instead of the github URL

```bash
jupyter-repo2docker --no-run \
    --image-name canine_bc \
    /Users/awatters/repos/FREYA
```

Note that since the notebooks are not in the mounted directory in the `docker run` above
any changes to the notebooks live within the container only and 
will not persist after the image is discarded.

In development mode to persistently edit the notebooks in place in your github repository 
working copy run the
container similar to this:

```bash
docker run -it --rm -p 8888:8888 \
   -v /Users/awatters/repos/FREYA/notebooks:/home/rstudio/notebooks:z \
   -v /Users/awatters/misc/kiley_graim/kgraim/data:/home/rstudio/notebooks/data:z \
   canine_bc:latest
```

The above command mounts the persistent directory containing scripts and notebooks
`/Users/awatters/repos/FREYA/notebooks` into
the container at the location
`/home/rstudio/notebooks`.  Changes made in the mounted directory will persist.

## Docker notes

If you do a lot of docker builds you will sometimes run out of space
in the docker build areas.  Use the docker
[`prune`](https://docs.docker.com/config/pruning/) command

```bash
docker image prune
```

to free up space or search for other similar methods.

## Todo (temp):

Fix notebook errors.

Output artifact (like diagrams) should go into a mounted persistent location.

Make repo public after vetting.

Provide mechanisms for uploading and downloading data that will work in Binder environment.
