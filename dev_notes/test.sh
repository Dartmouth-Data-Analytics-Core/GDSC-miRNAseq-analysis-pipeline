/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/bowtie2/bowtie2-2.4.4-linux-x86_64/bowtie2 \
    -x /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/human/hg38-bowtie2-index/hg38-bowtie2-index \
    -U seqcluster/A63_seqcluster.fastq.gz \
    -p 12 \
    --very-sensitive-local \
    --un seqcluster/unaligned_bowtie2.fastq.gz \
    -S seqcluster/seqclust_genome_aln.sam 2> seqcluster/seqclust_genome_aln.log

/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/samtools/samtools-1.11/samtools \
    view -Sb seqcluster/seqclust_genome_aln.sam | \
    /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/samtools/samtools-1.11/samtools sort -@ 4 - > seqcluster/seqclust_genome_aln.srt.bam 
/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/samtools/samtools-1.11/samtools \
    index seqcluster/seqclust_genome_aln.srt.bam


/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/shared-software/tools/bowtie/bowtie \
    -t -k 50 --best --strata -e 99999 --chunkmbs 2048 \
    --sam -p 8 -x /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/miRNA/miRNA_v2/dre/dre_hairpin_index/hairpins \
    --un seqcluster/hairpin_unaligned.fastq seqcluster/A63_seqcluster.fastq.gz \
    seqcluster/seqcluster_to_hairpin.aln.sam 2> seqcluster/seqcluster_to_hairpin.log

/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/samtools/samtools-1.11/samtools \
    view -Sb seqcluster/seqcluster_to_hairpin.aln.sam | \
    /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/samtools/samtools-1.11/samtools sort -@ 4 - > seqcluster/seqcluster_to_hairpin.aln.srt.bam 
/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/samtools/samtools-1.11/samtools \
    index seqcluster/seqcluster_to_hairpin.aln.srt.bam




#----- Build human bowtie1 index
mkdir hg38_bt1
bowtie-build GRCh38.primary_assembly.genome.fa hg38_bt1/genome


#----- Align to human genome
bowtie -a --best --strata -m 5000 \
  -x hg38_bt1/genome \
  seqcluster/D69_seqcluster.fastq.gz \
  --sam seqcluster/seqs.sam  
samtools view -Sb seqcluster/seqs.sam | samtools sort -@ 4- > seqcluster/seqs.srt.bam


# This is the exact mapping command used in nfcore smrnaseq
# but I think the seqcluster.fastq has been trimmed etc...
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

# When you pull the bam already made by nfsmrnaseq, THIS works....
# This does not work if we take the seqcluster fq and map to the human genome directly
mirtop gff \
    --add-extra \
    --sps hsa \
    --hairpin hsa_mirbase_hairpins.fa \
    --gtf hsa.gff3 \
    -o NEW_TEST \
    --format BAM \
    seqcluster/D69_test.bam # OUTPUT FROM NFSMRNASEQ

mirtop gff \
    --sps hsa \
    --hairpin /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/miRNA/miRNA_v2/hsa/hairpins_UtoT.fa \
    --gtf hsa.gff3 \
    -o NEW_TEST \
    seqcluster/seqs.srt.bam


#----- To clean the fasta file before indexing?????
bioawk \
        -c fastx '{gsub(/[^ATGCatgc]/, "N", $seq); sub(/ .*/, "", $name); print ">"$name"\n"$seq}' \
        genome.fa \
        > genome_clean.fa

# Reimplementation with regular awk
awk '
  BEGIN { seq=""; header="" }
  /^>/ {
    if (seq != "") {
      gsub(/[^ATGCatgc]/, "N", seq)
      print header
      print seq
    }
    header = $1       # take the first word of the header
    seq = ""
    next
  }
  {
    seq = seq $0
  }
  END {
    if (seq != "") {
      gsub(/[^ATGCatgc]/, "N", seq)
      print header
      print seq
    }
  }
' GRCh38.primary_assembly.genome.fa > genome_clean.fa


seqcluster \
    collapse \
    -m 1 --min_size 15 \
    -f D69.fastp.fastq.gz  \
    -o collapsed