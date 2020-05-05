
# Simulation Steps


### Simulate pedigrees

Run sim_code/simFamAJ.R with arguments 1-50 to generate families/simFamAJ_1.RData, ..., families/simFamAJ_50.RData 

$ sh sim_code/simFamAJ.sh
(1,000,000 families in total, 20,000 in each RData file)

To get perturbed versions of the families, run sim_code/simFamNoise.R with arguments 1-50 to generate families/simFamNoise_1.RData, ..., families/simFamNoise_50.RData 

$  sh sim_code/simFamAJ_noise.sh 



### Run BRCAPRO 

Run sim_code/run_brcapro_AJ.R with arguments 1-50 to generate brcapro_results/brcapro_AJ_1.RData, ..., brcapro_results/brcapro_AJ_50.RData. Run sim_code/run_brcapro_noise.R with arguments 1-50 to generate brcapro_results/brcapro_noise_1.RData, ..., brcapro_results/brcapro_noise_50.RData. 

$  cd sim_code 

$  sbatch --array=1-50 run_brcapro_AJ_array.sh 

$  sbatch --array=1-50 run_brcapro_noise_array.sh 

(run on the cluster)


### Standardize and flatten pedigrees

Run sim_code/nn_input.R to generate nn_input/nn_input.csv, nn_input/cnn_input.csv

Run sim_code/nn_input_noise.R to generate nn_input/nn_input_noise.csv, nn_input/cnn_input_noise.csv

Run sim_code/graph_mat.R to generate nn_input/graph_mat.csv (neighborhood matrix for CNN)


### Train NNs

Run the python scripts in training_code/fcnn_cnn.txt to generate results/*.csv (files with predictions prefixed by "pred", files with performance metrics prefixed by "res")


### Generate tables and figures

Run R scripts in tables_figures



# File Descriptions

### sim_code
- cgn.rel.counts.RData: CGN relative counts
- fam.gen.R: code for generating pedigrees with a user-specified structure
- geno.gen.R: code for generating BRCA1/BRCA2 genotypes for a pedigree
- pheno.gen.R: code for generating breast/ovarian cancer phenotypes for a pedigree with BRCA1/BRCA2 genotypes
- simFamAJ.R: code for generating 20,000 AJ families with BRCA1/BRCA2 genotypes and breast/ovarian cancer phenotypes (optional command line argument: seed value); each set of 20,000 families is saved to families/simFamAJ_*.RData
- simFamAJ.sh: generates 1,000,000 AJ families by running simFamAJ.R with seeds 1-50
- simFamAJ_noise.R: code for perturbing AJ families generated by simFamAJ.R via misreported diagnoses and ages (optional command line argument: seed value); each set of 20,000 perturbed families is saved to families_misrep/simFamNoise_*.RData
- simFamAJ_noise.sh: perturbs the 1,000,000 AJ families generated by simFamAJ.sh
- run_brcapro_AJ.R: runs BRCAPRO on families in families folder
- run_brcapro_AJ_array.sh: script for running BRCAPRO on correctly reported families cluster
- run_brcapro_noise.R: runs BRCAPRO on perturbed families in families_misrep folder
- run_brcapro_noise_array.sh: script for running BRCAPRO on perturbed families on cluster
- standardize_ped.R: code for standardizing pedigrees given a reference structure
- graph_mat.R: generates neighborhood matrix for CNN; saves matrix to nn_input/graph_mat.csv
- nn_input.R: standardizes families in families folder and generates NN inputs; NN inputs saved to nn_input/nn_input.RData, nn_input/nn_input.csv, nn_input/cnn_input.csv
- nn_input_noise.R: standardizes perturbed families in families_misrep folder and generates NN inputs; NN inputs saved to nn_input/nn_input_noise.RData, nn_input/nn_input_noise.csv, nn_input/cnn_input_noise.csv



### training_code
- graph_convolution.py: convolutional layer implementation from https://github.com/hechtlinger/graph_cnn/tree/master/code/graph_convolution.py
- graph_convolution_pro.py: implementation of CNN layer that passes on the proband's output to the next layer
- tune*.py: scripts for tuning FCNN and CNN hyperparameters; results saved in tuning_results
- tune*.sh: scripts for tuning FCNN and CNN hyperparameters; results saved in tuning_results
- plot*.py: scripts for plotting number of epochs vs validation loss / AUC / O/E (3 optional command line arguments: training set size, number of epochs, seed value); plots saved in tuning_results
- train_fcnn*.py: scripts for training FCNN (3 optional command line arguments: training set size, number of epochs, seed value); saves performance measures to results/res_*.csv and test set predictions to results/pred_*.csv 
- train_save_fcnn.py: script for training FCNN and saving the final model (3 optional command line arguments: training set size, number of epochs, seed value); saves model to training_code/fcnn*.h5, performance measures to results/res_*.csv, and test set predictions to results/pred_*.csv
- train_cnn*.py: scripts for training CNN (3 optional command line arguments: training set size, number of epochs, seed value); saves performance measures to results/res_cnn*.csv and test set predictions to results/pred_cnn*.csv 
- fcnn_cnn.txt: list of command line arguments for final NN models
- pred_fam.ipynb: evaluates test performance of NNs/LR/BRCAPRO across bootstrap replicates and calculates predictions for example family; performance measures from bootstrap replicates saved to results/auc_boot.csv, results/oe_boot.csv, results/bs_boot.csv, results/corr_boot.csv; predictions for example family saved to results/fam_pred.csv
- train_fcnn_noise.ipynb: evaluates test performance of NNs/LR/BRCAPRO across bootstrap replicates for perturbed families; performance measures from bootstrap replicates saved to results/auc_boot_noise2.csv, results/oe_boot_noise2.csv, results/bs_boot_noise2.csv


### tables_figures
- auc_plot.R: generate training set size vs AUC/correlation plot
- table_sim.R: generate performance table for correctly reported families
- table_sim_noise.R: generate performance table for misreported families
- family_interactions.R: generate NN inputs for example family; saved to nn_input/example_family*.csv
- fam_plot.R: generate plot of predictions for example family
- plot_ped.R: figure of example pedigree
- conv_nb.pptx: figure of CNN neighborhoods
- ped_map.tex: figure of pedigree standardization