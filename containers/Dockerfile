FROM continuumio/miniconda:4.6.14

# dependencies are:
#  samtools
#  python3
#  c++ boost libraries
#

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y

RUN conda install -c bioconda samtools
RUN pip install argparse
RUN pip install os.path
RUN pip install sys
RUN pip install subprocess

ADD ./wgbs_tools/bam2pat.py /usr/local/bin/
RUN chmod +x /usr/local/bin/bam2pat.py

