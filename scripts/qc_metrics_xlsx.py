import sys
import pandas as pd

#opening data
input_tsv = pd.read_csv(sys.argv[1], sep='\t')
#saving to xlsx
input_tsv.to_excel(sys.argv[2], index=False, sheet_name='Alignment_metrics')
