FROM python:3.10

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install libnode-dev librdf0 librdf0-dev unzip wget default-jre && apt-get clean  

RUN python -m pip install --upgrade pip

RUN wget https://github.com/boecker-lab/sirius/releases/download/v5.5.4/sirius-5.5.4-linux64-headless.zip && cd /usr/local/ && unzip /sirius-5.5.4-linux64-headless.zip && rm /sirius-5.5.4-linux64-headless.zip

ENV PATH="/usr/local/sirius/bin:${PATH}"

WORKDIR /opt/workdir
RUN mkdir /opt/workdir/data

COPY requirements.txt .
RUN pip3 install -r requirements.txt

COPY Workflow_Python_Functions.py .
COPY Workflow_Python_Script.py .