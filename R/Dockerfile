# To run the container one can use (assuming your data is in the current directory):
# docker run -v `pwd`:/opt/workdir/data zmahnoor/maw-r Rscript --no-save --no-restore --verbose Workflow_R_Script.r

FROM rocker/r-ver:4.2.0

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install libnode-dev librdf0 librdf0-dev unzip wget default-jre && apt-get clean

RUN apt-get install -y libnetcdf-dev libxml2-dev


RUN apt-get install -y \
  build-essential \
  libglpk40

RUN wget https://github.com/boecker-lab/sirius/releases/download/v4.9.12/sirius-4.9.12-linux64-headless.zip && cd /usr/local/ && unzip /sirius-4.9.12-linux64-headless.zip && rm /sirius-4.9.12-linux64-headless.zip

ENV PATH="/usr/local/sirius/bin:${PATH}"


WORKDIR /opt/workdir
RUN mkdir /opt/workdir/data

COPY install_packages.R .
RUN Rscript install_packages.R

COPY Workflow_R_Functions.r .
COPY Workflow_R_Script.r .
