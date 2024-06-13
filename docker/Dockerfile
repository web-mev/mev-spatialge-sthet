FROM ghcr.io/web-mev/base-spatialge-docker:sha-d38cc74eb7705991f074e32ce454e569965543ea

USER root
RUN mkdir -p /opt/resources && \
  chmod a+r /opt/resources && \
  chown $MAMBA_USER:$MAMBA_USER /opt/resources

USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1

ADD sthet.R /usr/local/bin
