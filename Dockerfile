FROM ghcr.io/prefix-dev/pixi:0.74.0 AS build

WORKDIR /opt/nf-rash

COPY pixi.toml pixi.lock ./
RUN pixi install --locked

FROM debian:bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/nrminor/nf-RASH"
LABEL org.opencontainers.image.title="nf-RASH"

COPY --from=build /opt/nf-rash/.pixi/envs/default /opt/nf-rash/.pixi/envs/default

ENV PATH="/opt/nf-rash/.pixi/envs/default/bin:${PATH}"

WORKDIR /work
