FROM python:3.9-bullseye

# Install conda
RUN apt-get update && \
    apt-get install -y build-essential gdal-bin libgdal-dev && \
    apt-get install -y wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Install conda env
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml

# # Install cron and start sheduling
# RUN apt-get update
# RUN apt-get -y install cron

# Copy code
WORKDIR /src
COPY . .

# Install package
RUN conda run -n dmpipeline python setup.py install

# Make sure the entry point script is executable
RUN chmod +x entrypoint.sh

# Run entrypoint when docker is started
ENTRYPOINT ["./entrypoint.sh"]