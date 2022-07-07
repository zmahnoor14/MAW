FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install libnode-dev librdf0 librdf0-dev python3-pip unzip wget default-jre && apt-get clean  

RUN apt-get remove r-base-core

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    rm -rf /var/lib/apt/lists/*

RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository ppa:marutter/rdev
RUN apt-get update
RUN apt-get upgrade
RUN apt-get install -y r-base
RUN apt-get install -y libnetcdf-dev


RUN wget https://github.com/boecker-lab/sirius/releases/download/v5.5.4/sirius-5.5.4-linux64-headless.zip && cd /usr/local/ && unzip /sirius-5.5.4-linux64-headless.zip && rm /sirius-5.5.4-linux64-headless.zip

ENV PATH="/usr/local/sirius/bin:${PATH}"

COPY requirements.txt .
RUN pip3 install -r requirements.txt

COPY install_packages.R .
RUN Rscript install_packages.R

WORKDIR /opt/workdir