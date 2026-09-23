bowtie \
    --threads 24 \
    --sam \
    -x /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/miRNA/miRNA_v2/hsa/new_hairpins/hairpins \
    -q \
    --un seqcluster/D69.unmapped.fastq \
    -t -k 50 --best --strata -e 99999 --chunkmbs 2048 \
    seqcluster/D69_seqcluster.fastq.gz \
    2>| >(tee D69.out >&2) \
    | samtools view  -@ 24 -bS -o seqcluster/D69_test.bam -

#Time loading forward index: 00:00:00
#Time loading mirror index: 00:00:00
#Seeded quality full-index search: 00:00:01
# reads processed: 6829
# reads with at least one alignment: 195 (2.86%)
# reads that failed to align: 6634 (97.14%)
#Reported 228 alignments
#Time searching: 00:00:01
#Overall time: 00:00:01


# When you pull the bam already made by nfsmrnaseq, THIS works....
# This does not work if we take the seqcluster fq and map to the human genome directly
mirtop gff \
    --add-extra \
    --sps hsa \
    --hairpin hsa_mirbase_hairpins.fa \
    --gtf hsa.gff3 \
    -o NEW_TEST \
    --format BAM \
    seqcluster/D69_test.bam

#/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/howe/miRNA-pipeline-v2-dev/bt1_env/.pixi/envs/default/lib/python3.12/site-packages/mirtop/libs/config.py:1: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
#  import pkg_resources
#/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/howe/miRNA-pipeline-v2-dev/bt1_env/.pixi/envs/default/lib/python3.12/site-packages/Bio/pairwise2.py:278: BiopythonDeprecationWarning: Bio.pairwise2 has been deprecated, and we intend to remove it in a future release of Biopython. As an alternative, please consider using Bio.Align.PairwiseAligner as a replacement, and contact the Biopython developers if you still need the Bio.pairwise2 module.
#  warnings.warn(
#['gff', '--add-extra', '--sps', 'hsa', '--hairpin', 'hsa_mirbase_hairpins.fa', '--gtf', 'hsa.gff3', '-o', 'NEW_TEST', '--format', 'BAM', 'seqcluster/D69_test.bam']
#02/27/2026 09:04:00 INFO Run annotation
#02/27/2026 09:04:00 INFO Reading seqcluster/D69_test.bam
#02/27/2026 09:04:01 INFO Hits: 6674
#02/27/2026 09:04:01 INFO Hits with indels 0
#02/27/2026 09:04:01 INFO Hits after clean: 6674
#02/27/2026 09:04:01 INFO Done.
#02/27/2026 09:04:01 INFO Valid hits (+/- reference miRNA): 0
#02/27/2026 09:04:01 INFO Skipped due to not precursor sequence: 0
#02/27/2026 09:04:01 INFO GFF miRNAs: 0
#02/27/2026 09:04:01 INFO GFF hits 0 by 0 reads
#02/27/2026 09:04:01 INFO Filtered by being duplicated: 0
#02/27/2026 09:04:01 INFO Filtered by being outside miRNA positions: 0
#02/27/2026 09:04:01 INFO Filtered by being low score: 0
#02/27/2026 09:04:01 INFO It took 0.005 minutes