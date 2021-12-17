# Imagene
ImaGene provides researchers with a transparent tool for which to begin radiogenomic analysis and explore possible further directions in their research. A web platform is available here: https://www.imagene.pgxguide.org/index.php

ImaGene comprises of four modules: a) data pre-processing, b) correlation Analysis, c) machine learning (ML), and d) reporting (Figure 1). The code for the software has been written in Python, and it utilizes several libraries such as scikit-learn (Pedregosa et al. 2011), matplotlib, seaborn and importr along with custom functions written to follow a systematic approach to analyze and pin point meaningful associations between imaging and omics features.

For analytical operations, ImaGene requires two input files, each in Comma-Separated Value (CSV) format – one containing imaging features (and their measurements) and another containing omics features for a set of tumor samples. The imaging features can be acquired from feature extraction software such as PyRadiomics, LIFEx, or RaCaT by processing the tumor images using the respective segmentation labels (Koçak et al. 2019; Pfaehler et al. 2019; van Griethuysen et al. 2017). The omics features can be acquired from studies conducted on tumor ROIs in biopsy samples processed in pathological laboratories, and may consist of data pertaining to gene expression, SV (including CNV), SNV, or DNA methylation scores.

REFERENCES:

Koçak, Burak, et al. (2019), 'Radiomics with artificial intelligence: A practical guide for beginners', Diagnostic and interventional radiology (Ankara, Turkey), 25 (6), 485-95.

Pedregosa, Fabian, et al. (2011), 'Scikit-learn: Machine learning in Python', Journal of machine learning research, 12, 2825-30.

Pfaehler, Elisabeth, et al. (2019), 'RaCaT: An open source and easy to use radiomics calculator tool', PLOS ONE, 14 (2), e0212223.

van Griethuysen, Joost J. M., et al. (2017), 'Computational Radiomics System to Decode the Radiographic Phenotype', Cancer research (Chicago, Ill.), 77 (21), E104-E07.




