import sys







#umi_dir='/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/snakemake/mirna-seq/DAC-miRNAseq-pipeline/git/DAC-miRNAseq-pipeline/umi_reads'
#mir_mapping_dir='/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/snakemake/mirna-seq/DAC-miRNAseq-pipeline/git/DAC-miRNAseq-pipeline/mirbase_alignment'
#genome_mapping_dir='/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/snakemake/mirna-seq/DAC-miRNAseq-pipeline/git/DAC-miRNAseq-pipeline/genome_alignment'

umi_dir = sys.argv[1]
mir_mapping_dir = sys.argv[2]
genome_mapping_dir = sys.argv[3]

import glob
umi_log_paths = glob.glob(umi_dir+'/*.log.txt')

#print (umi_log_paths)

umi_data_dict = {}
sample_list = []

for logfile in umi_log_paths:
    log = open(logfile,'r')

    sample_id = logfile.split('/')[-1].split('.')[0]
    sample_list.append(sample_id)

    for line in log:
        if line[0] == '#':
            continue

        sline = line.strip('\n')

        if 'INFO Input Reads:' in line:
            reads_before = int(sline.split(' ')[-1])
        if 'INFO Reads output:' in line:
            reads_after = int(sline.split(' ')[-1])

    reads_removed = reads_before - reads_after

#    print (reads_before, reads_after, reads_removed, sample_id)

    umi_data_dict[sample_id] = [reads_before, reads_removed, reads_after]

#print (umi_data_dict)


mir_mapping_log_paths = glob.glob(mir_mapping_dir+'/*mirbase.log.txt')
mirmap_data_dict = {}

for logfile in mir_mapping_log_paths:
    log = open(logfile,'r')

    sample_id = logfile.split('/')[-1].split('.')[0]

    for line in log:
        if 'reads; of these' in line:
            sline = line.strip()
            reads = int(sline.split(' ')[0])

        if 'aligned 0 times' in line:
            sline = line.strip()
            unaligned = int(sline.split(' ')[0])

    mirmap_data_dict[sample_id] = [reads-unaligned, unaligned]

#print (mirmap_data_dict)



genome_mapping_log_paths = glob.glob(genome_mapping_dir+'/*genome.log.txt')
genome_data_dict = {}

for logfile in genome_mapping_log_paths:
    log = open(logfile,'r')

    sample_id = logfile.split('/')[-1].split('.')[0].split('_')[:-1]
    sample_id = '_'.join(sample_id)

    for line in log:
        if 'reads; of these' in line:
            sline = line.strip()
            reads = int(sline.split(' ')[0])

        if 'aligned 0 times' in line:
            sline = line.strip()
            unaligned = int(sline.split(' ')[0])
            

    genome_data_dict[sample_id] = [reads-unaligned, unaligned]

#print (genome_data_dict)



#import pandas as pd


#df_umi = pd.DataFrame.from_dict(umi_data_dict)
#print(df_umi)

sample_list.sort()
print ('\t'.join(['Metric'] + sample_list))
print ('\t'.join(['# of reads'] + [str(umi_data_dict[x][0]) for x in sample_list]))
print ('\t'.join(['# of reads missing UMI'] + [str(umi_data_dict[x][1]) for x in sample_list]))
print ('\t'.join(['# of UMI-containing reads'] + [str(umi_data_dict[x][2]) for x in sample_list]))
print ('\t'.join(['# of reads mapping to mirna'] + [str(mirmap_data_dict[x][0]) for x in sample_list]))
print ('\t'.join(['# of reads mapping to genome'] + [str(genome_data_dict[x][0]) for x in sample_list]))
print ('\t'.join(['# of reads unaligned'] + [str(genome_data_dict[x][1]) for x in sample_list]))




