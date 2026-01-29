# Utilitiy Scripts
🛠️ Helpful day to day utilities that tie into GDSC-Pipelines and make life easier!


## Sample Sheet Generation
There are wo scripts that work in conjunction with one another:

1. `make_sample_sheet.sh`
2. `linkMeta.R`

The first script is a driver which generates a temporary sample sheet. The R-script makes use of `openxlsx` and `stringr` to open the slims metadata sheet from the sequencing facility and link file names to external IDs, preventing human error in making the sample sheet. By default, this pipeline outputs comma-separated files (.csv) but can be easily modified to output tab-separated files if desired in the `linkMeta.R` file. 

### Implementation
To view usage menu for the code, run the following:
```shell
bash make_sample_sheet.sh
```
The script takes two main arguments, the path to the raw data files on `GSR_Active` and a library layout (one of single or paired).
**By deafult, in the GSR folder, there should be an xlsx file named `metadata.xlsx`. Ensure this file is present before running the script!**

1. Run from **within** the `Utilities` folder by activating `sampleSheets` conda environment
```shell

cd Utilities

conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/sampleSheets

bash make_sample_sheet.sh /dartfs-hpc/rc/lab/G/GSR_Active/Labs/YourLab/YourProject paired

```

This will generate either `sample_fastq_list_single.csv` or `sample_fastq_list_paired.csv` **in the utilities folder** depending on your layout within your cloned DAC-RNASeq-Pipeline folder.
You can now specify this file name in `config.yaml` or in any of the configs in `prebuilt_configs`.

