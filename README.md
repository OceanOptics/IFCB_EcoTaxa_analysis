# IFCB_EcoTaxa_analysis
This code is a set of functions to process `.tsv` files exported from EcoTaxa for particle size and taxonomic group analysis. The functions were developed to analyze particle images collected using the Imaging FlowCytobot (IFCB; McLane Laboratories Inc.). The image files are uploaded to EcoTaxa (http://ecotaxa.obs-vlfr.fr) where the category of each image is predicted and may or may not be validated before the image information is exported in `.tsv` file format.

The code outputs a MATLAB structure which can then be used to plot the data and visualize particle size distributions and both the relative and absolute quantities of various particle and phytoplaknton categories. The plotting codes are not currently included in these functions.

## To get started
* Download the folder of functions, six total:
`run_extract_tsv_size_types.m`
`read_tsv_file.m`
`get_var_names.m`
`write_group_counts_file.m`
`bins_by_size.m`
`group_by_type.m`
  
* The driver code is `run_extract_tsv_size_types.m` and the complete input and output details are included in the header of the `.m` file, including two different example calls. The code is designed to provide the user with several options, inluding:
  * processing a single `.tsv` or a folder containing multiple `.tsv` files
  * providing a dilution factor to adjust the volume sampled, if necessary
  * choosing whether a `.txt` file containing the validated counts of each phytoplankton group is written or not
  * defining the bins for size analysis either by defining all bin edges (in microns), or by setting the min, max, and number of size bins resulting in a log-scale calculation of bin edges
  * using either volume-based or area-based diameter during size analyses
  
Please see the header of `run_extract_tsv_size_types.m` for further information on inputs and outputs. 
