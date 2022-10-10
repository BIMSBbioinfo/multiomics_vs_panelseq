# Multiomics vs Panelseq

Supporting repository for the manuscript of the same name. Here, one can find scripts and instructions to recreate the work partly or as a whole.


# Downloading the data
```get_data.sh``` script is here to access the data used in the work. Example usage: 
 ```bash
 git clone https://github.com/BIMSBbioinfo/multiomics_vs_panelseq.git
 cd multiomics_vs_panelseq
 bash get_data.sh -o <out_path>
 ```

# Building drug response prediction models

After the selected data were downloaded, one has multiple paths to follow. 
To generate models from the processed data the code section below would be sufficient:

 ```R
 usage:
 Rscript src/model_DR.R <dataset> <path_to_data_dir> <path_to_write_models> <number of runs> <number of cores>
 example (building models using CCLE or beatAML datasets, running 1 modeling run per drug using 10 cores:
 Rscript src/model_DR.R CCLE ./data ./results 1 10
 ```

 A modeling run for each drug can be repeated multiple times (we do this for PDX samples):
 ```R
 Rscript src/model_DR.R PDX ./data ./results 20 10
 ```

Once the models are built, performance metrics and variable importance metrics 
can be extracted and saved with the following
 ```R
 Rscript ./src/analyse_models.R ./results 20 #using 10 cores to process results
 ```

# Producing manuscript figures and tables 

Finally, figures and tables can be collated:
 ```R
 Rscript ./src/collate_figures_tables.R ./results
 ```


