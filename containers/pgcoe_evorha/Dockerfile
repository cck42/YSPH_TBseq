FROM continuumio/miniconda3:23.10.0-1

# chage to root dir, add conda to path
WORKDIR /root/
ENV PATH /opt/conda/envs/env/bin:$PATH

#install conda packages from environment (nb: linux64 arch)
COPY environment.yml .
RUN conda env update -n root -f environment.yml && conda clean -a

RUN conda install -y -c conda-forge openjdk

WORKDIR /root/
COPY evorha/evorha.jar /jars/evorha.jar

ENV PATH="${PATH}:/opt/miniconda/bin:/opt/miniconda/lib"
ENV LD_LIBRARY_PATH="${PATH}"
