FROM mambaorg/micromamba:1.3.1-alpine

USER root

RUN apk update && \
    apk upgrade && \
    apk add git && \
    apk add build-base && \
    apk add g++ && \
    apk add bsd-compat-headers && \
    apk add tiff

WORKDIR /app

COPY . /app

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG PIP_USE_PEP517=1

RUN micromamba install -y -n base -c conda-forge python=3.12
RUN pip install /app/.

RUN micromamba clean --all --yes

# Clean up the src code, as it is installed now
RUN rm -rf /app

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]