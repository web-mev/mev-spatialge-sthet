FROM ghcr.io/web-mev/base-spatialge-docker:sha-b9c9d80863af47a4ac85f420b96fd1ff57e2cde7

USER root
RUN mkdir -p /opt/resources && \
  chmod a+r /opt/resources && \
  chown $MAMBA_USER:$MAMBA_USER /opt/resources

USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# in addition to the packages in the base image, we need
# the annotation packages to enable creation of gene identifier mapping files
RUN micromamba install -y \
    -c conda-forge \
    -c bioconda \
    bioconductor-org.hs.eg.db=3.18.0 \
    bioconductor-org.mm.eg.db=3.18.0


#ADD create_gene_mappings.R /usr/local/bin

# Create the gene mapping files:
#RUN Rscript /usr/local/bin/create_gene_mappings.R human /opt/resources/human_genes.tsv
#RUN Rscript /usr/local/bin/create_gene_mappings.R mouse /opt/resources/mouse_genes.tsv

ADD sthet.R /usr/local/bin
ADD run.sh /usr/local/bin