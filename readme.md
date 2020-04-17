


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
