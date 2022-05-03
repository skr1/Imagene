# ImaGene
ImaGene provides researchers with a transparent tool for which to begin radiogenomic analysis and explore possible further directions in their research. A web platform is available here: https://www.imagene.pgxguide.org/index.php

ImaGene comprises of four modules: a) data pre-processing, b) correlation Analysis, c) machine learning (ML), and d) reporting (Figure 1). The code for the software has been written in Python, and it utilizes several libraries such as scikit-learn (Pedregosa et al. 2011), matplotlib, seaborn and importr along with custom functions written to follow a systematic approach to analyze and pin point meaningful associations between imaging and omics features.

For analytical operations, ImaGene requires two input files, each in Comma-Separated Value (CSV) format – one containing imaging features (and their measurements) and another containing omics features for a set of tumor samples. The imaging features can be acquired from feature extraction software such as PyRadiomics, LIFEx, or RaCaT by processing the tumor images using the respective segmentation labels (Koçak et al. 2019; Pfaehler et al. 2019; van Griethuysen et al. 2017). The omics features can be acquired from studies conducted on tumor ROIs in biopsy samples processed in pathological laboratories, and may consist of data pertaining to gene expression, SV (including CNV), SNV, or DNA methylation scores.

# New version 2:

Adding Imagene_v2_2.py (only runs from within the docker image which is packaged as "Imagene_v2_2.tar.gz" and attached above. Steps to load docker package and run it successfully aare described in the subsequent CLI-section below.

Changes in v2:
1. Automated permutations (n=20) conducted on labels that have AUC value > 0.9 and R-square > 0.25 from testing of models.
From the results of the permutation tests, users could calculate p_value and 95% Confidence Intervals for a) AUCs for labels at various decision     thresholds and b) R-squares for same labels.
   
2. Calculates R-squares using the sklearn.metrics.r2_score method. The results indeed match to the R-square calculated post ImaGene run using the formula presented in the manuscript i.e. R-square = 1 - v-square, where v-square = 1- ( RMSE(y)/Stdev(y) )^2, where 'RMSE/Stdev' is already indicated by ImaGene as one of its results for every label that is tested.

3. Calculating feature importances using "model.feature_importances_" method for DecisionTree model and "model.coef_" method for regression models. This indicates the importance (aka weightage) of each feature (i.e. imaging feature) when for example: a model is trained with imaging features to predict gene expressions of various genes (aka labels).

4. Added comments for function descriptions and other crucial steps (such as AUC calculations, feature importance reporting, etc.) within the respective functions.


# Command Line interface (CLI) operation using docker

Step 1: Install Docker on an Ubuntu Machine using the following docker-installation documentation: https://docs.docker.com/engine/install/ubuntu/
        For MacOS, use: https://docs.docker.com/desktop/mac/install/. Use docker using terminal app in Mac thereafter.

Step 2: Load the docker container package Imagene_v2_2.tar.gz
        
        docker load < Imagene_v2_2.tar.gz
        
Step 3: Check the loaded image

        docker images
        
        Output:
        --------------------------------------------------------------------------------------------------
        REPOSITORY              TAG                 IMAGE ID            CREATED             SIZE
        shreysukhadia/imagene   2.2                 37297fbf58cb        12 days ago         1.3 GB
        --------------------------------------------------------------------------------------------------

Step 4: Run the docker image

        docker run -v "$(pwd)":/data shreysukhadia/imagene:2.2 Imagene_v2_2.py --data /data/Supplementary_Table_BC_Radiomic_features.csv --label /data/Supplementary_Table_BC_Gene_FPKM.csv --config /data/config_IBC_LR.ini > log_file 2>&1 &
      
        Output: Same as indicated in supplementary files in the repo above and in the manuscript.

# In case user wants to try including patient outcomes alongwith radiomic and/or omic features:
The user could put the patient outcomes as a column in a comma-separated values (csv) file and use that as a “label” file, and put radiomics or genomics (or radiomics and genomics combined) features as columns in another csv file and use that as a “data” file. The rows in each of these files have to represent the sample names and the values in the feature columns would represent the measures of the respective feature (for example: texture, shape, size or intensity measure for radiomics, gene expressions for genomics and patient outcome decision measures for clinical outcomes). This way the AI models would get trained and tested to predict patient clinical outcomes from radiomics or genomics or radio-genomics combined features.



REFERENCES:

Koçak, Burak, et al. (2019), 'Radiomics with artificial intelligence: A practical guide for beginners', Diagnostic and interventional radiology (Ankara, Turkey), 25 (6), 485-95.

Pedregosa, Fabian, et al. (2011), 'Scikit-learn: Machine learning in Python', Journal of machine learning research, 12, 2825-30.

Pfaehler, Elisabeth, et al. (2019), 'RaCaT: An open source and easy to use radiomics calculator tool', PLOS ONE, 14 (2), e0212223.

van Griethuysen, Joost J. M., et al. (2017), 'Computational Radiomics System to Decode the Radiographic Phenotype', Cancer research (Chicago, Ill.), 77 (21), E104-E07.





