# Multi-stage build: rust+libclang for compilation, slim debian for runtime.
FROM rust:1.83-bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        clang \
        libclang-dev \
        cmake \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libssl-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
# Cache deps separately from source — re-copy of src/ only invalidates the
# dep-build layer when Cargo.toml/lock changes.
COPY Cargo.toml Cargo.lock* ./
RUN mkdir -p src && echo 'fn main() {}' > src/main.rs && \
    cargo build --release && \
    rm -rf src target/release/bam-readcount-rs target/release/deps/bam_readcount_rs*
COPY src ./src
RUN touch src/main.rs && cargo build --release && \
    strip target/release/bam-readcount-rs

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
        zlib1g \
        liblzma5 \
        libbz2-1.0 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/bam-readcount-rs /usr/local/bin/bam-readcount-rs

ENTRYPOINT ["/usr/local/bin/bam-readcount-rs"]
