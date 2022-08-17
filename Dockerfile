FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install libnode-dev librdf0 librdf0-dev python3-pip unzip wget default-jre && apt-get clean

RUN apt-get remove r-base-core

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    rm -rf /var/lib/apt/lists/*

RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN add-apt-repository ppa:marutter/rdev
RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y r-base
RUN apt-get install -y libnetcdf-dev


RUN apt-get install -y \
  build-essential \
  libglpk40

RUN wget https://github.com/boecker-lab/sirius/releases/download/v4.9.12/sirius-4.9.12-linux64-headless.zip && cd /usr/local/ && unzip /sirius-4.9.12-linux64-headless.zip && rm /sirius-4.9.12-linux64-headless.zip

ENV PATH="/usr/local/sirius/bin:${PATH}"


COPY requirements.txt .
RUN pip3 install -r requirements.txt

COPY install_packages.R .
RUN Rscript install_packages.R


WORKDIR /opt/workdir
