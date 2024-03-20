FROM continuumio/miniconda3:23.10.0-1

# chage to root dir, add conda to path
WORKDIR /root/
ENV PATH /opt/conda/envs/env/bin:$PATH

#install conda packages from environment (nb: linux64 arch)
COPY environment.yml .
RUN conda env update -n root -f environment.yml && conda clean -a

RUN conda install -y -c conda-forge openjdk
ENV PATH="${PATH}:/opt/miniconda/bin:/opt/miniconda/lib"
ENV LD_LIBRARY_PATH="${PATH}"


WORKDIR /root/
COPY evorha/evorha.jar /jars/evorha.jar

RUN conda install -c bioconda -c conda-forge bioconda::predicthaplo

# RUN conda install -y -c conda-forge make
# RUN conda install -y -c conda-forge gxx
# COPY PredictHaplo-2.1.4 /installs/PredictHaplo-2.1.4
# WORKDIR /installs/PredictHaplo-2.1.4
# RUN make 

# RUN conda install -y -c conda-forge -c defaults cmake
# COPY PredictHaplo-master /installs/PredictHaplo
# WORKDIR /installs/PredictHaplo
# RUN cmake -DCMAKE_BUILD_TYPE=Release -B build -S .
# RUN cmake --build build




