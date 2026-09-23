#! /usr/bin/env python

import sys
import numpy as np
import pandas as pd


def to_rpkm_pre_sliced(arr):
    tot_rc = np.sum(arr[:,1:], axis=0)
    return arr[:,1:] / np.outer(arr[:,0], tot_rc) * 10**9


def to_tpm_pre_sliced(arr):
    ind_read_per_base = arr[:,1:] / arr[:,0].reshape(arr[:,1:].shape[0], 1)
    return (10**6 * ind_read_per_base /
           np.sum(ind_read_per_base, axis=0).reshape(1, ind_read_per_base.shape[1]))


###########
# Extract data
###########

data = pd.read_csv(sys.argv[1], sep='\t')

path_name_rpkm = (sys.argv[1])[:-4] + "_rpkm.tsv" if sys.argv[2] == "single" else (sys.argv[1])[:-4] + "_fpkm.tsv"
path_name_tpm = (sys.argv[1])[:-4] + "_tpm.tsv"


#############################################################
# UPDATED COLUMN HANDLING
#############################################################

# Columns now:
# 0 Geneid
# 1 gene_name
# 2 gene_biotype
# 3 Chr
# 4 Start
# 5 End
# 6 Strand
# 7 Length
# 8+ counts

trimmed = np.array(data.iloc[:,7:].to_numpy(), dtype="float32")


###########
# Run normalization
###########

rpkm_np = to_rpkm_pre_sliced(trimmed)
tpm_np = to_tpm_pre_sliced(trimmed)


###########
# Convert back to pandas
###########

data_only_rpkm = pd.DataFrame(
    data=rpkm_np,
    index=data.index,
    columns=list(data.columns)[8:]
)

data_only_tpm = pd.DataFrame(
    data=tpm_np,
    index=data.index,
    columns=list(data.columns)[8:]
)


###########
# Reattach annotation columns
###########

rpkm_df = pd.concat([data.iloc[:,:8], data_only_rpkm], axis=1)
tpm_df = pd.concat([data.iloc[:,:8], data_only_tpm], axis=1)


###########
# Write output
###########

rpkm_df.to_csv(path_name_rpkm, sep="\t", index=False)
tpm_df.to_csv(path_name_tpm, sep="\t", index=False)

print("Finished RPKM/FPKM and TPM Normalization")