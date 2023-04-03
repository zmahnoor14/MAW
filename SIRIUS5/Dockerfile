# To run the container one can use (assuming your data is in the current directory):
# docker run -v `pwd`:/opt/workdir/data zmahnoor/maw-r Rscript --no-save --no-restore --verbose Workflow_R_Script.r

FROM rocker/r-ver:4.2.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install libnode-dev librdf0 librdf0-dev unzip wget default-jre && apt-get clean

RUN apt-get update

RUN apt-get install -y libnetcdf-dev libxml2-dev libcurl4-openssl-dev librsvg2-bin


RUN apt-get install -y \
  build-essential \
  libglpk40


WORKDIR /opt/workdir
RUN mkdir /opt/workdir/data


RUN wget https://github.com/boecker-lab/sirius/releases/download/v5.6.3/sirius-5.6.3-linux64-headless.zip

RUN unzip sirius-5.6.3-linux64-headless.zip -d /opt/workdir

ENV PATH="/opt/workdir/sirius/bin:${PATH}"

#COPY .env .

#RUN sirius login --user-env USER_ENV_NAME --password-env PASSWORD_ENV_NAME --show

COPY install_packages.R .
RUN Rscript install_packages.R


COPY Run_SIRIUS5.r .
COPY insilico/MS1DATA_SiriusP.tsv .

COPY insilico/SIRIUS/1_NA_iso_NA_MS1p_204.122756958008_SIRIUS_param.ms .
