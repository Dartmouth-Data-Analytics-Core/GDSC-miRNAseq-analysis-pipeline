#!/usr/bin/env bash
# Generate MD5 checksums for reference files defined in one or more pipeline
# config YAMLs. Reads reference paths directly from the config so nothing is
# hardcoded here.
#
# Usage:
#   bash scripts/generate_checksums.sh prebuilt_configs/human_config.yaml
#   bash scripts/generate_checksums.sh prebuilt_configs/human_config.yaml \
#                                       prebuilt_configs/mouse_config.yaml \
#                                       prebuilt_configs/zebrafish_config.yaml
#
# Output: prebuilt_configs/ref_md5s/{sps}_refs.md5 for each config.
# Validate: md5sum --check prebuilt_configs/ref_md5s/hsa_refs.md5

set -euo pipefail

BOWTIE1_EXTS=(.1.ebwt .2.ebwt .3.ebwt .4.ebwt .rev.1.ebwt .rev.2.ebwt)
BOWTIE2_EXTS=(.1.bt2  .2.bt2  .3.bt2  .4.bt2  .rev.1.bt2  .rev.2.bt2)

OUT_DIR="prebuilt_configs/ref_md5s"
mkdir -p "$OUT_DIR"

# --- helpers -----------------------------------------------------------------

yaml_val() {
    # Extract a scalar value from a YAML file: yaml_val key file
    # Returns empty string if key is absent or commented out.
    python3 -c "
import sys, yaml
data = yaml.safe_load(open('$2'))
print(data.get('$1', '') or '')
"
}

# --- process each config -----------------------------------------------------

for CONFIG in "$@"; do
    if [[ ! -f "$CONFIG" ]]; then
        echo "ERROR: config not found: $CONFIG" >&2
        exit 1
    fi

    SPS=$(yaml_val sps "$CONFIG")
    if [[ -z "$SPS" ]]; then
        echo "ERROR: 'sps' key missing from $CONFIG" >&2
        exit 1
    fi

    OUT="$OUT_DIR/${SPS}_refs.md5"
    echo "=== $CONFIG (sps=$SPS) ==="

    FILES=()
    MISSING=0

    # Single reference files
    for KEY in annotation_gtf hairpin_gff hairpin_fa spikein_reference_core; do
        VAL=$(yaml_val "$KEY" "$CONFIG")
        [[ -z "$VAL" ]] && continue
        FILES+=("$VAL")
    done

    # Bowtie2 index prefixes
    for KEY in padded_mature_index bowtie2_genome_index; do
        PREFIX=$(yaml_val "$KEY" "$CONFIG")
        [[ -z "$PREFIX" ]] && continue
        for EXT in "${BOWTIE2_EXTS[@]}"; do
            FILES+=("${PREFIX}${EXT}")
        done
    done

    # Bowtie1 index prefixes
    for KEY in bowtie1_hairpin_index bowtie1_mature_index bowtie_spikein_index; do
        PREFIX=$(yaml_val "$KEY" "$CONFIG")
        [[ -z "$PREFIX" ]] && continue
        for EXT in "${BOWTIE1_EXTS[@]}"; do
            FILES+=("${PREFIX}${EXT}")
        done
    done

    # Check all files exist before hashing
    for F in "${FILES[@]}"; do
        if [[ ! -f "$F" ]]; then
            echo "  MISSING: $F" >&2
            MISSING=1
        fi
    done
    if [[ $MISSING -eq 1 ]]; then
        echo "  Skipping $SPS — fix missing files above." >&2
        continue
    fi

    echo "  Hashing ${#FILES[@]} files..."
    md5sum "${FILES[@]}" > "$OUT"
    echo "  Written: $OUT"
    echo ""
done

echo "Done. To validate:"
for CONFIG in "$@"; do
    SPS=$(yaml_val sps "$CONFIG")
    echo "  md5sum --check $OUT_DIR/${SPS}_refs.md5"
done
