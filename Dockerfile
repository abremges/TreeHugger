FROM continuumio/miniconda3

ADD conda.yaml /tmp/conda.yaml
RUN conda env create -f /tmp/conda.yaml && \
    conda clean -ay && \
    echo "source activate treehugger" > ~/.bashrc
ENV PATH /opt/conda/envs/treehugger/bin:$PATH
