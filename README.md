# Multiomics vs Panelseq

Supporting repository for the manuscript of the same name. Here, one can find scripts and instructions to recreate the work partly or as a whole.

## Downloading the data
```get_data.sh``` script is here to access the data used in the work. Three download options are provided: ```raw, proc, all```. Example usage: 
 ```bash
 git clone https://github.com/BIMSBbioinfo/multiomics_vs_panelseq.git
 cd multiomics_vs_panelseq
 bash get_data.sh -s <raw|proc|all> -o <out_path>
 ```
## Preparing the data, building models and assessing predictions
After the selected data were downloaded, one has multiple paths to follow. 
To generate models from the processed data the code section below would be sufficient:
 ```R
 Rscript src/model_DR.R "CCLE" ./data <path_to_write_models>
 Rscript src/model_DR.R "PDX" ./data <path_to_write_models>
 ```
Alternatively, one can directly recreate the figures running this section:
 ```R
 Rscript src/plot_figures.R ./data <path_to_write_figures>
 ```