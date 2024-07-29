FROM ghcr.io/web-mev/base-spatialge-docker:sha-358218764bf7470b87f40095a1941c944ee5553f

USER root
RUN mkdir -p /opt/resources && \
  chmod a+r /opt/resources && \
  chown $MAMBA_USER:$MAMBA_USER /opt/resources

USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1

ADD sthet.R /usr/local/bin
