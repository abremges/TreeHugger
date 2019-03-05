FROM continuumio/miniconda3

ADD treehugger.yaml /tmp/treehugger.yaml
RUN conda env create -f /tmp/treehugger.yaml

RUN echo "source activate $(head -1 /tmp/treehugger.yaml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/treehugger.yaml | cut -d' ' -f2)/bin:$PATH
