# canine_breast_cancer_public
Public repository for canine breast cancer analysis

To build the docker image from the command line

```bash
jupyter-repo2docker --no-run \
    --image-name canine_bc \
    https://github.com/flatironinstitute/canine_breast_cancer_public
```

To run the image with a mounted data directory located at
`/Users/awatters/misc/kiley_graim/kgraim/data`

```bash
docker run -it --rm -p 8888:8888 \
   -v /Users/awatters/misc/kiley_graim/kgraim/data:/home/rstudio/work/data:z \
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
from the one shown above).

If you have trouble connecting to the jupyter server it may be because some other process
is responding to port 8888 on your system.  To find out whether processes are connected on
that port use

```bash
lsof -ti:8888 
```

Note that since the notebooks are not in the mounted directory in the `docker run` above
any changes to the notebooks will not persist after the image is discarded.

In development mode to persistently edit the notebooks in place in your github repository run the
container similar to this:

```bash
docker run -it --rm -p 8888:8888 \
   -v /Users/awatters/repos/canine_breast_cancer_public/work:/home/rstudio/work:z \
   -v /Users/awatters/misc/kiley_graim/kgraim/data:/home/rstudio/work/data:z \
   canine_bc:latest
```
