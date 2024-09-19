#! /usr/bin/env python

### Method for finding RPKM and TPM courtesy of:
### https://btep.ccr.cancer.gov/question/faq/what-is-the-difference-between-rpkm-fpkm-and-tpm/
### https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

####### 
## Method for TPM
#######

### For each individual read, let:
###     indiv_read_per_kilobase = indiv_read_count / (read_bases_count / 10**3)

### Then, we sum every individual read_per_kilobase value and divide by 10**6:
###     sum_read_per_base_rate_per_million = sum_over_reads{indiv_read_per_kilobase} / 10**6

### Finally, we divide an individual read's per kilobase value by the sum over one million:
###     indiv_transcript_per_million = indiv_read_per_kilobase / sum_read_per_base_rate_per_million

### It turns out that the 10**3's cancel out, so it is technically equivalent to:
###     indiv_transcript_per_million = (indiv_read_count / read_bases_count * 10**6) / (sum_over_reads{indiv_read_count / read_bases_count})

##### Input and Output
##### sys[1] is the tsv file name to take in and normalize, should be formatted in same way as featurecount's 
import sys
import numpy as np
import pandas as pd

# define TPM function 
def to_tpm(arr):
    trimmed = np.array(arr.iloc[:,1:].to_numpy(),dtype="float32")
    ind_read_per_base = trimmed[:,1:] / trimmed[:,0].reshape(trimmed.shape[0], 1)
    ind_read_per_base = np.nan_to_num(ind_read_per_base)
    return (10**6 * ind_read_per_base / np.sum(ind_read_per_base, axis = 0).reshape(1, ind_read_per_base.shape[1]))

########### Main script ###########

## Extract data
data = pd.read_csv(sys.argv[1], sep='\t')

# define output file path 
path_name_tpm = (sys.argv[1])[:-4]+"_tpm.tsv"

## Call TPM function
tpm_np = to_tpm(data)

## Convert back to Pandas for tsv file conversion
data_only_tpm = pd.DataFrame(data = tpm_np, index = data.index, columns = list(data.columns)[2:])

# add mirbaseID and length values back in 
tpm_df = pd.concat([data.iloc[:,:2], data_only_tpm], axis = 1)

# write to CSV 
tpm_df.to_csv(path_name_tpm, sep="\t", index=False)

print("Finished TPM Normalization")


rate=(data.iloc[:,2]/data.iloc[:,1])
rate = np.nan_to_num(rate)
((1548/23)/sum(rate))*10**6

(((1548/23)*10**3)/sum(rate))*10**6