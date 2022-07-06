FROM rocker/r-ver:4.2.0

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install libnode-dev librdf0 librdf0-dev unzip wget default-jre && apt-get clean

RUN apt-get install -y libnetcdf-dev

RUN wget https://github.com/boecker-lab/sirius/releases/download/v5.5.4/sirius-5.5.4-linux64-headless.zip && cd /usr/local/ && unzip /sirius-5.5.4-linux64-headless.zip && rm /sirius-5.5.4-linux64-headless.zip

ENV PATH="/usr/local/sirius/bin:${PATH}"

COPY install_packages.R .
RUN Rscript install_packages.R

WORKDIR /opt/workdir
