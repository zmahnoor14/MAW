FROM python:3.10

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install libnode-dev librdf0 librdf0-dev postgresql nvi unzip wget default-jre && apt-get clean  

RUN apt-get -y install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6

RUN python -m pip install --upgrade pip

RUN wget https://github.com/boecker-lab/sirius/releases/download/v5.5.4/sirius-5.5.4-linux64-headless.zip && cd /usr/local/ && unzip /sirius-5.5.4-linux64-headless.zip && rm /sirius-5.5.4-linux64-headless.zip

ENV PATH="/usr/local/sirius/bin:${PATH}"

WORKDIR /opt/workdir
RUN mkdir /opt/workdir/data

COPY requirements.txt .
RUN pip3 install -r requirements.txt

COPY basic_config.yml /etc


RUN sed -i -e 's/from collections import OrderedDict, Sequence/from collections import OrderedDict\nfrom collections.abc import Sequence/' /usr/local/lib/python3.10/site-packages/provenance/utils.py

COPY hmdb_dframe_str.csv /opt/workdir
ADD Example_Tyrosine /opt/workdir/Example_Tyrosine
ADD Example_MethionineSulfoxide /opt/workdir/Example_MethionineSulfoxide

COPY Workflow_Python_Functions.py .
COPY Workflow_Python_Script.py .

CMD sed -i -e "s/db\: postgresql\:\/\//db\: postgresql\:\/\/$POSTGRES_USER\:$POSTGRES_PASSWORD@database\/$POSTGRES_DB/" /etc/basic_config.yml

COPY requirements2.txt .
RUN pip3 install -r requirements2.txt
