FROM continuumio/miniconda3

ADD treehugger.yaml /tmp/treehugger.yaml
RUN conda env create -f /tmp/treehugger.yaml

RUN echo "source activate treehugger" > ~/.bashrc
ENV PATH /opt/conda/envs/treehugger/bin:$PATH
