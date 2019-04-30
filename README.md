# TCGA-TDA
Pipeline for the identification of novel driver genes using TDA on TCGA data used in:

Rabadan R, Mohamedi Y, Rubin U, Chu T, Elliott O, Ares L, Cal S, Obaya AJ, Levine AJ, and Camara PG, _"Identification of Relevant Genetic Alterations in Cancer using Topological Data Analysis"_. Submitted.

To run the pipeline follow the following steps:

1. Create a folder for the project (e.g. `LUAD`) and a subfolder where the RSEM expression files will be located (e.g. `LUAD/LUAD_RSEM`)
2. Download RSEM expression files and index file (`sdrf.txt`)
3. Run the script `tpm_matrix.R` with the following arguments:

```Rscript tpm_matrix.R -d "RSEM expression folder" -i "sdrf.txt file path" -p "project name" -a "Annotations.csv file path" -o "Anno_old_new.csv file path" -q "number of cores"```

An expression matrix will be generated after this step.

4. Download the corresponding Onconator file from Firehose Broad GDAC and run the script `maf_process3.R` with the following arguments:

```Rscript maf_process3.R -f "path to Onconator tar.gz file" -p "project name"```

Mutations matrices will be generated after this step.

5. Run the script `Big_matrix.R` with the following arguments:

```Rscript Big_matrix.R -e "path to expression matrix FULL_TPM_MATRIX" -b "path to mutation matrix FULL_MUTATIONS_BINARY" -s "path to mutation matrix FULL_MUTATIONS_SYNONYMOUS" -n "path to mutation matrix FULL_MUTATIONS_NON_SYNONYMOUS" -p "project name"```

.csv and .h5 files will be generated after this step.

6. Run the python script `ayasdi_test.py` with the following arguments:

```python ayasdi_test.py "project name" "path to .csv file generated in step 5"```

7. To assess the mutational load run the script `connectivity8.R` with the following arguments:

```Rscript connectivity8.R -p "number of permutations" --mutload TRUE -m "path to .h5 file generated in step 5"```

8. Run the script `connectivity8.R` with the following arguments:

```Rscript connectivity8.R -p "number of permutations" -m "path to .h5 file generated in step 5" -z FALSE -q "number of cores" -t "mutation frequency threshold" -r "mutation downsampling threshold in log10 scale" --maf "path to .maf file" -g 350```
