#getting base image 
FROM danielsoler93/pele:v1.0.2
MAINTAINER Daniel Soler <daniel.soler@nostrumbiodiscovery.com>
#Set env variables for PELE
ENV PELE="/PELE-repo/"
ENV SCHRODINGER="/opt/"
ENV PATH="/usr/local/sbin:/usr/local/bin:/usr/local/bin/miniconda3/bin/:/usr/sbin:/usr/bin:/sbin:/bin"
#Install dependencies
RUN conda install -c conda-forge -c nostrumbiodiscovery pele_platform==1.4.0 --yes
#Create user
#RUN groupadd -g 999 appuser && \
#    useradd -r -u 999 -g appuser appuser
RUN useradd -ms /bin/bash appuser
USER appuser
WORKDIR /home/appuser
#Define build arguments
ARG PDB=""
ARG INPUT=""
# Copy input files
RUN echo $PDB
COPY $PDB "/home/appuser/"
COPY $INPUT "/home/appuser/"
#Launch job
RUN python -m pele_platform.main $INPUT
