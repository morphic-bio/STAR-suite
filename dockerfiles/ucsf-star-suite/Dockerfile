FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        coreutils \
        findutils \
        gzip \
        libgomp1 \
        libssl3 \
        pigz \
        python3 \
        seqtk \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/star-suite
COPY core/legacy/source/STAR.release /opt/star-suite/core/legacy/source/STAR.release
COPY scripts/paper/run_ucsf_corrected_production_workflow.sh /opt/star-suite/scripts/paper/run_ucsf_corrected_production_workflow.sh
COPY scripts/run_ucsf_perturb_yremove_batch.sh /opt/star-suite/scripts/run_ucsf_perturb_yremove_batch.sh
COPY scripts/lib/star_trim_qc.sh /opt/star-suite/scripts/lib/star_trim_qc.sh
COPY scripts/ucsf_parity/render_star_inputs_from_cr_config.py /opt/star-suite/scripts/ucsf_parity/render_star_inputs_from_cr_config.py
COPY scripts/prepare_velocyto_mex.py /opt/star-suite/scripts/prepare_velocyto_mex.py

RUN chmod +x \
      /opt/star-suite/core/legacy/source/STAR.release \
      /opt/star-suite/scripts/paper/run_ucsf_corrected_production_workflow.sh \
      /opt/star-suite/scripts/run_ucsf_perturb_yremove_batch.sh \
      /opt/star-suite/scripts/ucsf_parity/render_star_inputs_from_cr_config.py \
      /opt/star-suite/scripts/prepare_velocyto_mex.py

ENV STAR_BIN=/opt/star-suite/core/legacy/source/STAR.release
ENTRYPOINT ["/bin/bash"]
