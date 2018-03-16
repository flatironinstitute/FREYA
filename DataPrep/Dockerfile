FROM broadinstitute/java-baseimage

RUN apt-get update \
    && apt-get install -y bzip2 git libbz2-dev liblzma-dev make parallel python-scipy python-pysam python-setuptools samtools vim zip zlib1g-dev

# HTSeq will be installed into the system python area
WORKDIR /tmp
RUN wget https://github.com/simon-anders/htseq/archive/release_0.9.1.zip \
    && unzip -a release_0.9.1.zip \
    && rm release_0.9.1.zip
WORKDIR /tmp/htseq-release_0.9.1
RUN python setup.py install
WORKDIR /tmp
RUN rm -r htseq-release_0.9.1

# Layout of workflow components (there are other files mixed in with these):
# /usr/cmwf/
#  DEXSeq/inst/python_scripts/dexseq_count.py
#  disBatch/disBatch.py
#  FastQC/fastqc
#  hisat2-2.0.4/hisat2
#  VCF_mutation_picker.0.5.py
# /usr/cmwf/jars/
#  GenomeAnalysisTK-3.8-0.jar
#  picard.jar
#  snpEff/snpEff.jar
RUN mkdir -p /usr/cmwf/jars
WORKDIR /usr/cmwf/jars

# GATK
RUN wget -O gatk.tar.bz2 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-0-ge9d806836' \
    && tar jxf gatk.tar.bz2 \
    && mv GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar GenomeAnalysisTK-3.8-0.jar \
    && rm -r gatk.tar.bz2 GenomeAnalysisTK-3.8-0-ge9d806836

# picard
RUN wget https://github.com/broadinstitute/picard/releases/download/2.17.4/picard.jar

# snpEff
RUN wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3p_core.zip && unzip -a snpEff_v4_3p_core.zip && rm snpEff_v4_3p_core.zip

WORKDIR /usr/cmwf

# dexseq_count.py
RUN wget http://bioconductor.org/packages/3.4/bioc/src/contrib/DEXSeq_1.20.2.tar.gz \
    && tar zxf DEXSeq_1.20.2.tar.gz \
    && rm DEXSeq_1.20.2.tar.gz

# disBatch.py
RUN git clone --recursive https://github.com/flatironinstitute/disBatch.git
RUN disBatch/disBatch.py --fix-paths

# The next two blocks derive in part from NBISweden/K9-WGS-Pipeline files.
# fastqc
ENV FASTQC_VERSION=0.11.5
RUN wget --quiet -O fastqc_v${FASTQC_VERSION}.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip \
    && unzip fastqc_v${FASTQC_VERSION}.zip \
    && chmod ugo+x FastQC/fastqc \
    && rm fastqc_v${FASTQC_VERSION}.zip

# hisat2
ENV HISAT2_URL=ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-source.zip
ENV HISAT2_ARCHIVE_FILE hisat.zip
ENV HISAT2_BUILD_DIR hisat2-2.0.4
RUN wget -qO $HISAT2_ARCHIVE_FILE $HISAT2_URL \
    && unzip $HISAT2_ARCHIVE_FILE \
    && rm $HISAT2_ARCHIVE_FILE \
    && cd  $HISAT2_BUILD_DIR \
    && make -j 4 

# VCF_mutation_picker. TODO: Is "raw" access to github fragile?
RUN wget https://raw.github.com/flatironinstitute/FREYA/master/DataPrep/VCF_mutation_picker.0.6.py
