FROM rocker/r-ver:4.2.0

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install libnode-dev librdf0 librdf0-dev unzip wget default-jre && apt-get clean

<<<<<<< HEAD
RUN apt-get install -y libnetcdf-dev libxml2-dev
=======
RUN apt-get install -y libnetcdf-dev
>>>>>>> refs/remotes/origin/main


RUN apt-get install -y \
  build-essential \
  libglpk40

RUN wget https://github.com/boecker-lab/sirius/releases/download/v4.9.12/sirius-4.9.12-linux64-headless.zip && cd /usr/local/ && unzip /sirius-4.9.12-linux64-headless.zip && rm /sirius-4.9.12-linux64-headless.zip

ENV PATH="/usr/local/sirius/bin:${PATH}"

COPY install_packages.R .
RUN Rscript install_packages.R

COPY Workflow_R_Functions.r .

WORKDIR /opt/workdir