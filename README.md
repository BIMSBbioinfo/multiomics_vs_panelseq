# Multiomics vs Panelseq

Supporting repository for the manuscript of the same name. Here, one can find scripts and instructions to recreate the work partly or as a whole.

## Analysing CCLE and PDX datasets 

###  Downloading the data
```get_data.sh``` script is here to access the data used in the work. Three download options are provided: ```raw, proc, all```. Example usage: 
 ```bash
 git clone https://github.com/BIMSBbioinfo/multiomics_vs_panelseq.git
 cd multiomics_vs_panelseq
 bash get_data.sh -s <raw|proc|all> -o <out_path>
 ```
### Preparing the data, building models and assessing predictions
After the selected data were downloaded, one has multiple paths to follow. 
To generate models from the processed data the code section below would be sufficient:
 ```R
 Rscript src/model_DR.R "CCLE" <path_to_data_dir> <path_to_write_models>
 Rscript src/model_DR.R "PDX" <path_to_data_dir> <path_to_write_models>
 ```
Alternatively, one can directly recreate the figures running this section:
 ```R
 Rscript src/plot_figures.R <path_to_data_dir> <path_to_write_figures>
 ```


## Analysing BeatAML dataset

Below is a description of how to download, process and analyse BeatAML dataset from the [orcestra.ca database](https://www.orcestra.ca/pset/canonical).
The outputs of the scripts can be found under ./data. 

	- data/beatAML.prepared.RDS (output from step 2)
	- data/beatAML.stats.tsv and data/beatAML.plot.pdf (output from step 3)

1. First create a working folder:

```
 mkdir BeatAML 
 cd BeatAML
```

2. Preparing the data (assuming the current folder is BeatAML)
  
This will import BeatAML_2018 dataset from [orcestra.ca database](https://www.orcestra.ca/pset/canonical) and process the 
gene expression, mutation, and drug sensitivity data to make it ready for downstream analysis.  
  
  ```
  Rscript ../src/beatAML_analysis/prepare_data.R ../data 
  ```
The prepared dataset can be found in ./data/beatAML.prepared.RDS. 

3. Modeling drug-sensitivity using genomic/transcriptomic features 

This will build two sets of models for every drug that was treated against at least 100 samples. 
One model using only mutation features for genes in the onkokb panel, and a second model 
using gene-set scores derived from rna-seq dataset along on top of the mutation features from the first model. 
A train/test split will be created for each drug (70/30); and a repated-cross-validation procedure will 
be followed on the training data. The performance of the models will be tested on the test split and 
statistics will be generated for each drug. 

The generated stats are saved in beatAML.stats.tsv, and a summary figure is created in beatAML.plot.pdf. 

  ```
  Rscript ../src/beatAML_analysis/analysis.R
  ```







