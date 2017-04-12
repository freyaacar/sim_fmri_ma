# sim_fmri_ma
The purpose of this project is to simulate an ecologically valid fmri meta-analysis.  

To be able to test meta-analytical tools and methods it is important to have a meta-analysis set that resembles a real meta-analysis in within- and between "study" variance but in the meantime does not suffer from the same lack of data as fMRI meta-analyses. In this project we aim to develop such a meta-analysis set and test it's within- and between "study" variance and the influence of family ties on the variability.  

We use a meta-analysis on pain which can be downloaded from neurovault (http://neurovault.org/collections/1425/) (https://www.nature.com/articles/sdata2016102). If the folder with the data is put into the working directory with the code, everything should work fine.  

## Done
* Read in t-maps  
* Compute effect sizes  

## To do
* Calculate within- and between study variance  
* Determine other parameters (number of peaks, clusters, peak heights, cluster sizes, sample size, ...)  
* Find new datasets (+ extract as many parameters as possible)
* Export different formats (t-maps, ES maps, SPM output, FSL output, peak coordinates)
* Compare these parameters to the parameters of a meta-analysis constructed from HCP dataset
* Compute between-study variance of HCP dataset
* Add more real meta-analyses/studies
