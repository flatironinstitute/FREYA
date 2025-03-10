FROM rocker/tidyverse:3.4.2

RUN apt-get -y remove bash-completion

RUN apt-get update && \
    apt-get -y install python3-pip bash-completion && \
    pip3 install --no-cache-dir notebook==5.4.1 && \
    apt-get purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

#RUN apt-get -y install -f bash-completion

ENV NB_USER rstudio
ENV NB_UID 1000
ENV HOME /home/rstudio
WORKDIR ${HOME}

USER ${NB_USER}

# Set up R Kernel for Jupyter
RUN R --quiet -e "install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest'))"
RUN R --quiet -e "devtools::install_github('IRkernel/IRkernel')"
RUN R --quiet -e "IRkernel::installspec()"

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
RUN R -f /home/rstudio/packages.R
USER root
RUN chown -R ${NB_UID}:${NB_UID} ${HOME}
RUN chsh -s /bin/bash ${NB_USER}
RUN sudo apt-get install python3-pip
RUN python3 -m pip install jp_proxy_widget
RUN python3 -m pip install requests
RUN jupyter nbextension enable --py --sys-prefix jp_proxy_widget

USER ${NB_USER}

# Run install.r if it exists
RUN if [ -f install.r ]; then R --quiet -f install.r; fi

CMD jupyter notebook --ip 0.0.0.0
