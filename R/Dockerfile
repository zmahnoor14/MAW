# To run the container one can use (assuming your data is in the current directory):
# docker run -v `pwd`:/opt/workdir/data zmahnoor/maw-r Rscript --no-save --no-restore --verbose Workflow_R_Script.r

FROM rocker/r-ver:4.2.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install \
	libxml2-dev \
	librsvg2-bin \
	libssl-dev \
	libnode-dev \
	librdf0 \
	librdf0-dev \
	libnetcdf-dev \
	libharfbuzz-dev \
	libfribidi-dev \
	libgit2-dev \
	libfreetype6-dev \
	libpng-dev \
	libtiff5-dev \
	libjpeg-dev \
	libfontconfig1-dev \
	python3.10 \
	python3-pip \
	unzip \
	wget \
	default-jre \
	build-essential \
	libglpk40 \
	ant \
	default-jdk \
	git \
	&& apt-get clean

<<<<<<< HEAD
RUN apt-get update

RUN apt-get install -y libnetcdf-dev libxml2-dev libcurl4-openssl-dev librsvg2-bin libssl-dev

RUN apt-get update

RUN apt-get -y install libcurl4-gnutls-dev 

RUN apt-get update

RUN apt-get install -y \
  build-essential \
  libglpk40
=======
RUN pip3 install provbook
>>>>>>> provenance

RUN apt-get install -y libharfbuzz-dev libfribidi-dev libgit2-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libfontconfig1-dev

RUN apt-get install -y python3.10
RUN apt-get install -y python3-pip

RUN pip3 install jupyter
RUN pip3 install provbook


RUN wget https://github.com/boecker-lab/sirius/releases/download/v4.9.12/sirius-4.9.12-linux64-headless.zip && cd /usr/local/ && unzip /sirius-4.9.12-linux64-headless.zip && rm /sirius-4.9.12-linux64-headless.zip

ENV PATH="/usr/local/sirius/bin:${PATH}"

WORKDIR /opt/workdir
RUN mkdir /opt/workdir/data

COPY install_packages.R .
RUN Rscript install_packages.R

COPY install_packages_devtools.R .
RUN Rscript install_packages_devtools.R

<<<<<<< HEAD
RUN apt-get update

RUN apt-get install -y ant default-jdk git

RUN cd /opt/workdir

RUN git clone https://github.com/End-to-end-provenance/DDG-Explorer.git

WORKDIR /opt/workdir/DDG-Explorer

RUN ant build-project

RUN ant ddg-explorer

COPY gnps.rda /opt/workdir/
COPY mbankNIST.rda /opt/workdir/
COPY hmdb.rda /opt/workdir/
COPY Example_Tyrosine.mzML /opt/workdir/
COPY Example_MethionineSulfoxide.mzML /opt/workdir/
WORKDIR /opt/workdir
COPY Workflow_R_Functions.r .
COPY Workflow_R_Script.r .
COPY Workflow_R_Script.ipynb .
=======
RUN git clone https://github.com/End-to-end-provenance/DDG-Explorer.git
WORKDIR /opt/workdir/DDG-Explorer
RUN ant build-project
RUN ant ddg-explorer

>>>>>>> provenance
