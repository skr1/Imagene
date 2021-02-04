# Imagene
A multitude of laboratories across the globe conduct radiogenomic studies employing a variety of statistics and machine learning methods internally without revealing the code and parameters used for such analysis, limiting reproducibility [1]and robustness of such studies.   At present, the field lacks a comprehensive and robust software able to intake both imaging and omics features of a tumor as input and perform statistical correlation and machine learning studies to unveil associations between them and cast them into learned models [2, 3]. Moreover, there is a need for an algorithmic approach with a flexibility in parameter setting for to model associations between these features to derive robust conclusions. To conceal this gap, we have developed Imagene, a software that is first in its kind that adopts an algorithmic approach to integrate statistics with AI for the robust radiogenomic analysis of tumors in patients. Imagene facilitates researchers with a systematic radiogenomic analysis  with flexible parameters and a comprehensive HTML report depicting the analytical steps enhanced with the corresponding metrics and visualization items. Overall, Imagene promises to improve reproducibility and transparency, and robustness of the radiogenomics analysis through an easy to operate command-line (CLI) python-coded software, a flexible parameter configuration file, and an intuitive and easy-to-use report that could be shared between the research groups. Imagene could be further improved by the continuous contribution from scientific community.

Imagene constitutes four modules: a) Data pre-processing and normalization, b) Statistical correlations, c) Machine learning (ML) and d) Report generation. It requires user to prepare and input three files, first being a user config file with all parameter settings (Table 1), second a comma-separated values (csv) formatted imaging-feature file with samples names as rows and features as columns. The values for these features may exist on varying numerical scales that are handled by data preprocessing and normalization step based on the normalization parameter mentioned in the config file. Thirdly, a similar Omics feature containing csv file. Imaging features could be pre-acquired by processing tumor scans with the respective segmentation labels through any of the feature extraction software as listed in the introduction section above. Omics features could be any including gene expressions, CNVs, SNVs, DNA methylation scores, etc. Imagene is coded completely in Python language, uses R, Matlab and scikit-learn libraries, and could be operated easily mentioning appropriate arguments on the command line.

    Configuration Parameters:	Values
    Correlation Method	default: Spearman (options: Pearson)
    Correlation Threshold	default: 0.5 (flexible)
    pValue correction method	default:BH (options: holm, hochberg, hommel, bonferroni, BH, BY, fdr)
    Data Type	default:Imaging (option: Gene)
    Label Type	Gene (option: Imaging)
    Train size	0.9 (flexible)
    Test size	0.1 (flexible)
    Normalization method	default: Stand_scaler (options: min_max , Stand_scaler , zscore, MaxAbsScaler)
    Mode	default:Train (options: Train, validate, predict)
    Model Type	default: Decision Tree (options: Linear Regression, Linear Model, LASSO)
    Cross Validation splitter	cv: 2 (flexible)
    <Section for other model parameters for each model type>	Could be filled with as many parameters possible or kept empty.
    
References:

1. Brito, L. (2020) Recommendations to enhance rigor and reproducibility in biomedical research. Gigascience. [Online] 9 (6), .
2. Colen, F. (2014) NCI Workshop Report: Clinical and Computational Requirements for Correlating Imaging Phenotypes with Genomics Signatures. Translational oncology. [Online] 7 (5), 556–569.
3. Bodalal, T. (2019) Radiogenomics: bridging imaging and genomics. Abdominal radiology (New York). [Online] 44 (6), 1960–1984.

