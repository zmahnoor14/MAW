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

RUN pip3 install provbook

RUN wget https://github.com/boecker-lab/sirius/releases/download/v4.9.12/sirius-4.9.12-linux64-headless.zip && cd /usr/local/ && unzip /sirius-4.9.12-linux64-headless.zip && rm /sirius-4.9.12-linux64-headless.zip

ENV PATH="/usr/local/sirius/bin:${PATH}"

WORKDIR /opt/workdir
RUN mkdir /opt/workdir/data

COPY install_packages.R .
RUN Rscript install_packages.R

COPY install_packages_devtools.R .
RUN Rscript install_packages_devtools.R

RUN git clone https://github.com/End-to-end-provenance/DDG-Explorer.git
WORKDIR /opt/workdir/DDG-Explorer
RUN ant build-project
RUN ant ddg-explorer

