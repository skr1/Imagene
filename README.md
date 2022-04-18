# ImaGene
ImaGene provides researchers with a transparent tool for which to begin radiogenomic analysis and explore possible further directions in their research. A web platform is available here: https://www.imagene.pgxguide.org/index.php

ImaGene comprises of four modules: a) data pre-processing, b) correlation Analysis, c) machine learning (ML), and d) reporting (Figure 1). The code for the software has been written in Python, and it utilizes several libraries such as scikit-learn (Pedregosa et al. 2011), matplotlib, seaborn and importr along with custom functions written to follow a systematic approach to analyze and pin point meaningful associations between imaging and omics features.

For analytical operations, ImaGene requires two input files, each in Comma-Separated Value (CSV) format – one containing imaging features (and their measurements) and another containing omics features for a set of tumor samples. The imaging features can be acquired from feature extraction software such as PyRadiomics, LIFEx, or RaCaT by processing the tumor images using the respective segmentation labels (Koçak et al. 2019; Pfaehler et al. 2019; van Griethuysen et al. 2017). The omics features can be acquired from studies conducted on tumor ROIs in biopsy samples processed in pathological laboratories, and may consist of data pertaining to gene expression, SV (including CNV), SNV, or DNA methylation scores.


# Command Line interface (CLI) using docker

Step 1: Install Docker on an Ubuntu Machine using the following docker-installation documentation: https://docs.docker.com/engine/install/ubuntu/
        For MacOS, use: https://docs.docker.com/desktop/mac/install/. Use docker using terminal app in Mac thereafter.

Step 2: Load the docker container package Imagene_v2.tar.gz
        
        docker load < Imagene_v2.tar.gz
        
Step 3: Check the loaded image

        REPOSITORY              TAG                 IMAGE ID            CREATED             SIZE
        shreysukhadia/imagene   2.0                 f9cbd4bab876        2 days ago          1.3 GB
        
Step 4: Run the docker image

        docker run -v "$(pwd)":/data shreysukhadia/imagene:2.0 Imagene_v2.py --data /data/Supplementary_Table_BC_Radiomic_features.csv --label /data/Supplementary_Table_BC_Gene_FPKM.csv --config /data/config_IBC_LR.ini > log_file 2>&1 &
      


REFERENCES:

Koçak, Burak, et al. (2019), 'Radiomics with artificial intelligence: A practical guide for beginners', Diagnostic and interventional radiology (Ankara, Turkey), 25 (6), 485-95.

Pedregosa, Fabian, et al. (2011), 'Scikit-learn: Machine learning in Python', Journal of machine learning research, 12, 2825-30.

Pfaehler, Elisabeth, et al. (2019), 'RaCaT: An open source and easy to use radiomics calculator tool', PLOS ONE, 14 (2), e0212223.

van Griethuysen, Joost J. M., et al. (2017), 'Computational Radiomics System to Decode the Radiographic Phenotype', Cancer research (Chicago, Ill.), 77 (21), E104-E07.






