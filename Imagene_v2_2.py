##Imagene.py
##Note: This script is only operable when within the docker container package. Check github for more info on running the package.
##Version 1.02: ##Added if condition to check if the low ratio RMSE plot file exists before writing it to HTML file
##Version 1.03.1: ## a) Changed the conditions for binarization of feature columns in Y_test and Y_pred for AUC calculations. Now, if value<threshold then it gets binarized to 0 else to 1.
##              ## b) Changed decision threshold list to [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0].
##              ## c) Edited mp.xlim and ylim parameters for AUC v/s decision_threshold plot to include left=0.0 and right=1.0 values (for xlim).
##Version 1.04: ## Changed the decision threshold list to [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9].
##              ## Changed the if condition in 1038 to match length of the AUC_value_list and decision_threshold_list to >=8. Made the correspoding change in print statement in line 1046 accordingly.
##              ## Changed xlim: now left=0.1 and right=0.9
##              ## Reverted the condition for binarization of feature columns in Y_test and Y_pred for AUC calculations. Now, if value<=threshold then it gets binarized to 0 else to 1.
##Version 1.05: ##Using pickle module instead of joblib module to load model.pkl
##Version 1.1: ##Adding multiTask lasso and multitask elastic net models
##Version 2.0: Added Feature importances and permutation tests for R-square and AUC. R-square gets calculated using sklearn.metrics.r2_score method.
##Version 2.2: if pVal_adj_method="none" (i.e. p_adjust=p_value) or corr_threshold (aka correlation coefficient threshold) < 0.0 (for ex: -1.0), then no <0.05 filtering on p_adjust value.
##Note: Please monitor github commits for version updates after Version: 2.2.
##Author: Shrey Sukhadia
#!/usr/bin/python
from cProfile import label
import matplotlib as mpl
mpl.use('Agg')
import os, re, sys, math
import argparse
import numpy as np
import pandas as pd
import joblib
import pickle
import configparser
import ast
import argparse
from sklearn.feature_selection import SelectFromModel
from sklearn import preprocessing
from sklearn.preprocessing import MaxAbsScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import Normalizer
from sklearn.preprocessing import binarize
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import auc
from scipy import stats
from scipy.interpolate import make_interp_spline, BSpline
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn import linear_model
from sklearn import tree
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import MultiTaskElasticNet
from sklearn.linear_model import LinearRegression
from sklearn import metrics
from sklearn.metrics import make_scorer
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import seaborn as sb
import matplotlib.pyplot as mp
import scikitplot as skplt
from math import sqrt
import base64
from datetime import datetime
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

##CORRELATION function establishes correlations between gene and imaging features and calculates FDR adjusted pvalues further to establish statistical significance of those correlations.
def correlation(d_1,d_2,d_1_header,d_2_header, model_type, corr_method, corr_threshold,pVal_adjust_method, tagDir):
    '''
        Conduct correlations for the given data and label features depending on the correlation method (aka pearson or spearman) and correlation co-efficient threshold chosen by user.
        Returns signficantly(p_adjust<0.05) correlated features. Users have option to select the method used to adjust p-values. Check parameter settings below.
        ##Note that all parameters are set in respective config.ini (if running imagene through CLI on linux/mac terminal) or directly on ImaGene's web platform if running there.

        :param d1: data features
        :type d1: dataframe
        :param d2: label features
        :type d1: dataframe
        :param d2: label features
        :type d1: dataframe
        :param d_1_header: data header
        :type d_1_header: list
        :param d_2_header: label header
        :type d_2_header: list
        :param model_type: Type of the model, options: {'DecisionTree','LinearRegression', 'LinearModel' , 'LASSO', 'multiTaskLASSO', 'multiTaskLinearModel'}, default: None 
                            ##Note: This parameter is mainly used here to tag output file name with appropriate model_type, i.e. for which ML-build execute later in BuildModel function.
                            ##Helps to track experimental output files better, i.e. when correlated features (output from this function) are fed further to building the ML models in the BuildModel call from Process function defined towards the bottom of the script.
        :type model_type: str
        :param corr_method: Correlation method used, options: {'Pearson', 'Spearman'}
        :type corr_method: str
        :param corr_threshold: Threshold for correlation co-efficient, options: {-1.0,0.0,0.1...1.0}, default:0.5
                            ##Note: To silent correlation threshold filtering, set the threshold to -1.0
        :type corr_threshold: float
        :param pVal_adjust_method: The method to correct p-values, options: {'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'}, default: BH
        :type pVal_adjust_method: str
        :param tagDir: The sub-directory name to store files in. Could be kept empty (i.e. "") or filled with a name. Gets filled automatically when pipeline runs in second mode (automatically), i.e. "FeatureSelection".
                        ##Note: This is not a user-governed parameter. Internally set and used by the pipeline.
        :type tagDir: str
    '''
    pval=dict()
    pcorr=dict()
    Var1_list=[]; Var2_list=[]; spCorr_list=[]; pvalue_list=[]
    for i in d_1_header:
        for j in d_2_header:
            if(corr_method=="spearman"):
                pCorr,pv=stats.spearmanr(d_1[i], d_2[j])
            elif(corr_method=="pearson"):
                pCorr,pv=stats.pearsonr(d_1[i], d_2[j])
            else:
                pCorr,pv=stats.spearmanr(d_1[i], d_2[j])
            Var1_list.append(i); Var2_list.append(j); spCorr_list.append(pCorr); pvalue_list.append(pv)
    rstats = importr('stats')

    pAJ_list = rstats.p_adjust(FloatVector(pvalue_list), method = pVal_adjust_method)
    
    
    Master_df = pd.DataFrame({'Var1': Var1_list, 'Var2': Var2_list, 'correlation': spCorr_list, 'p_value': pvalue_list, 'p_adjust': pAJ_list })
    
    
    Master_df_sorted=Master_df.sort_values(by=['p_adjust'])
    Master_df_sorted.to_csv("/data/All_correlations.csv")

    if(pVal_adjust_method=="none" or corr_threshold < 0.0):
        Master_df_sorted_sgn=Master_df_sorted
    else:
        Master_df_sorted_sgn=Master_df_sorted[Master_df_sorted['p_adjust']<0.05]
    Master_df_sorted_sgn_filtered=Master_df_sorted_sgn[abs(Master_df_sorted_sgn['correlation'])>corr_threshold]
    
    
    
    Master_df_sorted_sgn.to_csv("/data/Significant_correlations.csv")
    Master_df_sorted_sgn_filtered.to_csv("/data/Significant_correlations_gt_corr_threshold.csv")
    List_of_Var1=Master_df_sorted_sgn_filtered["Var1"].tolist()
    List_of_Var2=Master_df_sorted_sgn_filtered["Var2"].tolist()
    fC_List_Var1=[]
    for i in List_of_Var1:
        fC=i.split("_")[0]
        fC_List_Var1.append(fC)
    fC_uniq = sorted(set(fC_List_Var1))

    image_tag_list=[]

    for i in fC_uniq:
        print(i)
        fC_df=Master_df_sorted_sgn_filtered[Master_df_sorted_sgn_filtered['Var1'].str.contains(i, regex=False, na=False)]
        fC_df_pivot=fC_df.pivot_table(index="Var1",columns="Var2",values="correlation",fill_value=0)
        print("For "+i+" features:")
        print(fC_df_pivot)
        try:
            sb.clustermap(fC_df_pivot)
        except ValueError as err:
            mp.clf()
            sb.heatmap(fC_df_pivot)
        
        if(pVal_adjust_method=="none" or corr_threshold < 0.0):
            mp.title("Correlations for "+i+" features", size=16)
        else:
            mp.title("Top significant correlations (FDR_adjusted_pValue<0.05) for "+i+" features", size=16)
        mp.xticks(rotation=90)
        mp.yticks(size=7)
        mp.savefig('/data/'+'Correlation_for_'+i+'_features.png',orientation='landscape',dpi=90,bbox_inches='tight')
        data_image = open("/data/"+tagDir+'Correlation_for_'+i+'_features.png', 'rb').read().encode('base64').replace('\n', '')
        image_tag_list.append('<img src="data:image/png;base64,{0}" style="max-width:50%;">'.format(data_image))
        mp.clf()

    outfileHTML=open("/data/"+tagDir+model_type+".output.html",'a')
    outfileHTML.write("<h2 style=text-align:center;>"+"--------------------------Multivariate Correlations ("+corr_method+" based)-----------------------"+"</h2>"+"\n".join(image_tag_list)+"\n")
    outfileHTML.close()

    return(sorted(set(List_of_Var1)),sorted(set(List_of_Var2)))

##READ datasets
def read_dataset(dataset):
    '''
        Reads dataset and converts to dataframe.
        Returns dataframe.

        :param dataset: dataset file (".csv")
        :type dataset: file
    '''
    if os.path.isfile(dataset):
        raw_dataframe = pd.read_csv(dataset, sep=',')
    else:
        print("{} should be csv/tsv file ".format(dataset))
        sys.exit()
    dataframe = raw_dataframe.copy()
    print("Shape of dataframe:{}".format(dataframe.shape))
    return dataframe

##NORMALIZATION function
def normal_dataframe(dataframe, norm_type, header, normalize = True):
    '''
        Normalizes data based on normalization method specified by user.
        Returns normalized dataframe.

        :param dataframe: dataframe for the data to be normalized
        :type dataframe: dataframe
        :param norm_type: Normalization type/method, options: {'min_max' , 'Stand_scaler' , 'MaxAbsScaler', 'none'}, default: None
        :type norm_type: str
        :param header: A header-list containing column names for dataframe.
        :type dataframe: list
        :param normalize: A boolean parameter indicating whether to normalize the dataframe or not. default: True.
        :type normalize: bool
    '''
    if normalize:
        if norm_type =='min_max':
            scaler = MinMaxScaler()
            dataframe_scaled = scaler.fit_transform(dataframe)
        elif norm_type == 'Stand_scaler':
            scaler = StandardScaler()
            dataframe_scaled = scaler.fit_transform(dataframe)
        elif norm_type == 'MaxAbsScaler':
            scaler = MaxAbsScaler()
            dataframe_scaled = scaler.fit_transform(dataframe)
        else:
            print("Invalid normalization type/method detected. Skipping normalization. If you wish to normalize this dataset, then correct the normalization_method in the config file and rerun. Proceeding with no normalization")
            return dataframe
        dataframe_scaled=pd.DataFrame(data=dataframe_scaled,columns=header)
        return dataframe_scaled
    else:
        return dataframe

##PREPROCESSING of datasets    
def preprocessing(dataframe , label, data_type, label_type, mode, tagDir, checkNA = True):
    '''
        Preprocessing of dataframe to check SampleIDs and eliminate columns containing 'nan' value in one or more than one cells.
        Returns processed dataframe.

        :param dataframe: dataframe for the data to be preprocessed.
        :type dataframe: dataframe
        :param label: dataframe for the label(s) to be preprocessed.
        :type label: dataframe
        :param data_type: A string indicating the type of data, ex: "imaging" or "gene".
        :type data_type: str
        :param label_type: A string indicating the type of label, ex: "imaging" or "gene".
        :type label_type: str
        :param mode: Mode of run, options: {'Train','predict','validate'}. ##Note: This parameter is no more in this function and would be removed from argument list in next version.
        :type mode: str
        :param tagDir: The sub-directory name to store files in. Could be kept empty (i.e. "") or filled with a name. Gets filled automatically when pipeline runs in second mode (automatically), i.e. "FeatureSelection".
                        ##Note: This is not a user-governed parameter. Internally set and used by the pipeline.
        :type tagDir: str
        :param checkNA: A boolean parameter indicating whether to check for nan values in dataframe. default: True
        :type checkNA: bool
    '''
    outfileHTML=open("/data/"+tagDir+model_type+".output.html",'a')
    outfileHTML.write("<h3>"+"No. of "+data_type+" features provided: "+str(len(dataframe.columns)-1)+"</h3>")
    if(isinstance(label,pd.DataFrame)):
        outfileHTML.write("<h3>"+"No. of "+label_type+" features provided:"+str(len(label.columns)-1)+"</h3>")   
    if checkNA:
        dataframe.replace("", np.nan, inplace=True)
        if dataframe.isnull().values.any():
            dataframe = dataframe.dropna(axis=1, how='any')
        if isinstance(label,pd.DataFrame):
            label.replace("", np.nan, inplace=True)
            if label.isnull().values.any():
                label = label.dropna(axis=1, how='any')
                #print label
    if isinstance(label,pd.DataFrame):
        if dataframe['ID'].equals(label['ID']):
            outfileHTML.write("<h3>"+"SampleID check results: 'The SampleIDs match for "+data_type+" and "+label_type+" features'"+"</h3>"+"\n")
            #outfileHTML.write("<h3>"+"No. of samples: "+"</h3>"+"\n")
            print("The SampleIDs match for "+data_type+" and "+label_type+" features")
            #dataframe = dataframe.drop(['ID'],axis = 1)
            #label = label.drop(['ID'],axis = 1)
        else:
            sys.exit("The SampleIDs in data and label vary. Cannot proceed further. Please fix and rerun!")
        label = label.drop(['ID'],axis = 1)
        label_header=list(label.keys())
    else:
        label_header="NA"
    sampleIDs = dataframe['ID']
    outfileHTML.write("<h3>"+"No. of samples: "+str(len(sampleIDs))+"</h3>"+"\n")
    dataframe = dataframe.drop(['ID'],axis = 1)
    dataframe_header=list(dataframe.keys())
    
    outfileHTML.close()

    print(dataframe)
    print(dataframe.shape[1])
    
    if isinstance(label,pd.DataFrame):
        print(label)
        print(label.shape[1])
    
    return dataframe , label, sampleIDs, label_header, dataframe_header

##SPLITTING of data into TRAIN and TEST datasets   
def splitdata(dataframe , label, t_size, mode_, data_normalize_method, label_normalize_method, data_type, label_type, dataframe_header, label_header, tagDir):
    '''
        Split data (and labels) into train and test sets and calls normalization for them.
        Returns normalized train and test data and labels.

        :param dataframe: dataframe for the data.
        :type dataframe: dataframe
        :param label: dataframe for the label(s).
        :type label: dataframe
        :param data_type: A string indicating the type of data, ex: "imaging" or "gene".
        :type data_type: str
        :param label_type: A string indicating the type of label, ex: "imaging" or "gene".
        :type label_type: str
        :param mode: Mode of run, options: {'Train','predict','validate'}. ##Note: This parameter is no more in this function and would be removed from argument list in next version.
        :type mode: str
        :param tagDir: The sub-directory name to store files in. Could be kept empty (i.e. "") or filled with a name. Gets filled automatically when pipeline runs in second mode (automatically), i.e. "FeatureSelection".
                        ##Note: This is not a user-governed parameter. Internally set and used by the pipeline.
        :type tagDir: str
        :param data_normalize_method: Method used to normalize data, options: {'min_max' , 'Stand_scaler' , 'MaxAbsScaler', 'none'}, default: None
        :type data_normalize_method: str
        :param label_normalize_method: Method used to normalize labels, options: {'min_max' , 'Stand_scaler' , 'MaxAbsScaler', 'none'}, default: None
        :type label_normalize_method: str
        :param data_type: A string indicating the type of data, ex: "imaging" or "gene".
        :type data_type: str
        :param label_type: A string indicating the type of label, ex: "imaging" or "gene".
        :type label_type: str
        :param data_header: data header
        :type data_header: list
        :param label_header: label header
        :type label_header: list
    '''
    outfileHTML=open("/data/"+tagDir+model_type+".output.html",'a')
    train, test , Y_train , Y_test = train_test_split(dataframe, label , test_size = t_size)
    
    ##Converting numpy arrays to dataframe to perform normalization on them
    train=pd.DataFrame(data=train,columns=dataframe_header)
    print(train)
    test=pd.DataFrame(data=test,columns=dataframe_header)
    print(test)
    Y_train=pd.DataFrame(data=Y_train,columns=label_header)
    print(Y_train)
    Y_test=pd.DataFrame(data=Y_test,columns=label_header)
    print(Y_test)

    ##Introducing normalizations for TRAIN and TEST datasets.
    if data_normalize_method != 'none':
        print("performing "+data_normalize_method+" normalization for "+data_type+" features")
        
        ##NORMALIZING TRAIN data
        outfileHTML.write("<h3>"+"performing "+data_normalize_method+" normalization for "+data_type+" features for TRAIN set"+"</h3>"+"\n")
        train= normal_dataframe(train, data_normalize_method, dataframe_header)
        print("Printing Normalized Train data:")
        print(train)

        ##NORMALIZING TEST set
        outfileHTML.write("<h3>"+"performing "+data_normalize_method+" normalization for "+data_type+" features for TEST set"+"</h3>"+"\n")
        test= normal_dataframe(test, data_normalize_method, dataframe_header)
        print("Printing Normalized Test data:")
        print(test)

    if isinstance(label,pd.DataFrame):
        if label_normalize_method != 'none':
            print("performing "+label_normalize_method+" for "+label_type+" features")
            
            ##NORMALIZING TRAINING label
            outfileHTML.write("<h3>"+"performing "+label_normalize_method+" for "+label_type+" features for TRAIN set"+"</h3>"+"\n")
            Y_train = normal_dataframe(Y_train, label_normalize_method, label_header)
            print("Printing Normalized Train label:")
            print(Y_train)

            ##NORMALIZING TEST label
            outfileHTML.write("<h3>"+"performing "+label_normalize_method+" for "+label_type+" features for TEST set"+"</h3>"+"\n")
            Y_test = normal_dataframe(Y_test, label_normalize_method, label_header)
            print("Printing Normalized Test label:")
            print(Y_test)
    
    outfileHTML.write("<h1 style=text-align:center;color:purple>"+"-------------------------------Number of Samples for Training and Testing---------------------------------"+"</h1>"+"\n")
    outfileHTML.write("<h3>"+"No. of samples for training:{}".format(len(train))+"</h3>"+"\n")
    outfileHTML.write("<h3>"+"No. of samples for test:{}".format(len(test))+"</h3>"+"\n")
    outfileHTML.close()
    print("Trainig data:{} , Testing data:{} ".format(len(train) ,len(test)))

    return train , Y_train , test , Y_test

##MODEL BUILDING: multiple options of modeltypes accepted from user
def BuildModel(train , Y_train , test , Y_test , method, params, cv_par, scoring_par, gridsearch, param_grid, select_label_var_list, select_data_var_list, data_type, label_type, featureSelFrmModel_flag, tagDir, trainmodel):
    '''
        Initializing the model and executing its training. Calls evaluate method to test the test dataset.
        Returns the model

        :param train: dataframe for train data.
        :type train: dataframe
        :param Y_train: dataframe for train labels.
        :type Y_train: dataframe
        :param test: dataframe for test data.
        :type test: dataframe
        :param Y_test: dataframe for test labels.
        :type Y_train: dataframe
        :param method: Model type, options: {'DecisionTree','LinearRegression', 'LinearModel' , 'LASSO', 'multiTaskLASSO', 'multiTaskLinearModel'}, default: None
        :type method: str
        :param params: Model parameters. If empty, default parameters used per scikit-learn library.
        :type params: dict
        :param cv_par: K-fold Cross validation splitter number, default: 2
        :type cv_par: int
        :param scoring_par: Score metric, default: neg_mean_square_error
        :type scoring_par: str
        :param gridsearch: Indicating whether to perform gridsearch for training the model, options: 'True' or 'False', default: 'False'
        :type gridsearch: str
        :param param_grid: Parameters for gridsearch for the given model type.
        :type param_grid: dict
        :param data_type: A string indicating the type of data, ex: "imaging" or "gene".
        :type data_type: str
        :param label_type: A string indicating the type of label, ex: "imaging" or "gene".
        :type label_type: str
        :param featureSelFrmModel_flag: flag to run 'FeatureSelection' mode.
        :type featureSelFrmModel_flag: bool
        :param tagDir: The sub-directory name to store files in. Could be kept empty (i.e. "") or filled with a name. Gets filled automatically when pipeline runs in second mode (automatically), i.e. "FeatureSelection".
                        ##Note: This is not a user-governed parameter. Internally set and used by the pipeline.
        :type tagDir: str
        :param select_data_var_list: List of data features to train.
        :type select_data_var_list: list
        :param select_label_var_list: List of label features to train.
        :type select_label_var_list: list
        :param data_type: A string indicating the type of data, ex: "imaging" or "gene".
        :type data_type: str
        :param label_type: A string indicating the type of label, ex: "imaging" or "gene".
        :type label_type: str
        :param trainmodel: Whether to train a model or not, options: True, False, default: None
        :type trainmodel: bool
    '''
    outfileHTML=open("/data/"+tagDir+model_type+".output.html",'a')
    outfileHTML.write("<h2 style=text-align:center;color:blue>"+"--------------------------Model Summary-----------------------"+"</h2>"+"\n")
    if method in ['DecisionTree','LinearRegression', 'LinearModel' , 'LASSO', 'multiTaskLASSO', 'multiTaskLinearModel']:
        if gridsearch == 'True':
            if 'cv' in param_grid.keys():
                cv_grid=param_grid['cv']
                del param_grid['cv']
            else:
                cv_grid=None
            if 'scoring' in param_grid.keys():
                scoring_grid=param_grid['scoring']
                del param_grid['scoring']
            else:
                scoring_grid=None
        if method == 'LASSO':
            if gridsearch == 'True':
                try:
                    print(" Starting grid search for LASSO")
                    model = GridSearchCV(linear_model.Lasso(), param_grid=param_grid,cv=cv_grid,scoring=scoring_grid)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                model = linear_model.Lasso(**params)
        elif method == 'DecisionTree':
            if gridsearch == 'True':
                try:
                    print("starting grid search for Decision Tress")
                    model = GridSearchCV(tree.DecisionTreeRegressor(), param_grid=param_grid, scoring=scoring_grid, cv=cv_grid)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                model = tree.DecisionTreeRegressor(**params)  
        elif method == 'LinearRegression':
            if gridsearch == 'True':
                try:
                    print("starting grid search for Linear Regression")
                    model = GridSearchCV(LinearRegression(), param_grid=param_grid, scoring=scoring_grid, cv=cv_grid)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                model = LinearRegression(**params)
        elif method == 'LinearModel': 
            if gridsearch == 'True':
                try:
                    print("Starting grid search for ElasticNet")
                    model = GridSearchCV(ElasticNet(), param_grid=param_grid , cv = cv_grid , scoring=scoring_grid)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                model = ElasticNet(**params)
        elif method == 'multiTaskLASSO':
            if gridsearch == 'True':
                try:
                    print(" Starting grid search for MultiTaskLASSO")
                    model = GridSearchCV(linear_model.MultiTaskLasso(), param_grid=param_grid,cv=cv_grid,scoring=scoring_grid)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                model = linear_model.MultiTaskLasso(**params)
        elif method == 'multiTaskLinearModel': 
            if gridsearch == 'True':
                try:
                    print("Starting grid search for MultiTaskElasticNet")
                    model = GridSearchCV(MultiTaskElasticNet(), param_grid=param_grid , cv = cv_grid , scoring=scoring_grid)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                model = MultiTaskElasticNet(**params)
    else:
        print("options are :DecisionTree, LinearRegression, LinearModel, LASSO, LinearModel (aka ElasticNet), multiTaskLinearModel, multiTaskLASSO")
    if trainmodel:
        '''
        performing training
        '''
        ## ADDING SelectFromModel for Feature Selection using model and using only those features for training and testing further (change featureSelFrmModel_flag to sfm_flag later)
        if(featureSelFrmModel_flag==1):
            column_headers__=train.columns
            selector = SelectFromModel(estimator=model).fit(train, Y_train)
            feature_selected_or_not_=selector.get_support()
            ##Make a dataframe of that array which has the header as names of each feature.
            fsd=pd.DataFrame(data=feature_selected_or_not_)
            fsd=pd.DataFrame.transpose(fsd)
            fsd.columns=column_headers__
            print(fsd)
            ##Drop the features that are "False", i.e. not selected.
            fsd=fsd.drop(columns=fsd.columns[(fsd == False).any()])
            print("These are the features selected by SelectFromModel function")
            print(fsd)
            fsd.to_csv("/data/"+tagDir+'_'+model_type+'_features_selected.txt')
            feature_headers__=fsd.columns
            ##Converting numpy array to dataframe with headers for selected features.
            train=train[feature_headers__]
            print("This is the train set post feature selection")
            print(train)
            ##Selecting same features in test as well.
            print("This is the test set post feature selection")
            test=test[feature_headers__]
            print(test)
        outfileHTML.write("<h3>"+"Model Type : "+method+"</h3>"+"\n")
        
        ##EXECUTING GRID-SEARCH based on whether the grid_search parameter set True
        if gridsearch == 'True':
            grid_result = model.fit(train ,Y_train)
            outfileHTML.write("<h4>"+"Grid Search Metrics"+"</h4>"+"\n")
            outfileHTML.write("\n"+"<h5>"+'Best Score : '+str(grid_result.best_score_)+"</h5>"+"\n")
        else:
            scores = cross_val_score(model.fit(train ,Y_train), train , Y_train, cv=cv_par , scoring = scoring_par)
            outfileHTML.write("<h4>"+"Cross Validation Metrics:"+"</h4>"+"\n")
            outfileHTML.write("<h5>"+"Parameters: cv="+str(cv_par)+"  scoring="+str(scoring_par)+"</h5>"+"\n")
            outfileHTML.write("<h5>"+"Cross validation score:{}".format(-1*scores.mean())+"</h5>"+"\n")

        store_params=model.get_params();
        outfileHTML.write("<h3>"+"Model Parameters:"+"</h3>"+"\n")
        for i in store_params.keys():
            outfileHTML.write("<h4>"+str(i)+":"+str(store_params[i])+"</h4>")
        outfileHTML.close()
        
        ##CALL EVALUATION for TRAIN set
        evaluate(model,train,Y_train,select_label_var_list,'train_eval',model_type, data_type, label_type, tagDir, gridsearch)
        
        print(Y_test)
        ##CALL EVALUATION for TEST set
        returned_label_resulted=evaluate(model,test,Y_test,select_label_var_list,'test_eval',model_type, data_type, label_type, tagDir, gridsearch)
        
        ##PERMUT each label column of the test dataset 20 times and validate the model. This aids in p-value calculation post run.
        for l in returned_label_resulted:
            for n in range(1,21):
                Y_test[[l]]=np.random.permutation(Y_test[[l]])
                evaluate(model, test, Y_test , select_label_var_list, 'validation'+'_permut_'+str(n)+'_'+str(l), model_type, data_type, label_type, tagDir, gridsearch)
        
    return model

##EVALUATE the models for TEST set or validation set.
def evaluate(model , test , Y_test, select_label_var_list, prefix, model_type, data_type, label_type, tagDir, gridsearch):
    '''
        Evaluating the model preformance and reporting variety of metrics
        Returns labels that got AUC>0.9 and R-square>0.25

        :param test: data for model evaluation
        :type test: dataframe
        :param Y_test: labels for model evaluation 
        :type Y_test: dataframe
        :param model: Model to evaluate
        :type model: model
        :param model_type: Model type, options: {'DecisionTree','LinearRegression', 'LinearModel' , 'LASSO', 'multiTaskLASSO', 'multiTaskLinearModel'}, default: None
        :type model_type: str
        :param gridsearch: Indicating whether to perform gridsearch for training the model, options: 'True' or 'False', default: 'False'
        :type gridsearch: str
        :param prefix: Prefix for output file-name
        :type prefix: str
        :param data_type: A string indicating the type of data, ex: "imaging" or "gene".
        :type data_type: str
        :param label_type: A string indicating the type of label, ex: "imaging" or "gene".
        :type label_type: str
        :param select_label_var_list: List of label features
        :type select_label_var_list: list
        :param tagDir: The sub-directory name to store files in. Could be kept empty (i.e. "") or filled with a name. Gets filled automatically when pipeline runs in second mode (automatically), i.e. "FeatureSelection".
                        ##Note: This is not a user-governed parameter. Internally set and used by the pipeline.
        :type tagDir: str
    '''
    outfileHTML=open("/data/"+tagDir+model_type+".output.html",'a')
    Y_pred = model.predict(test)
    print(Y_test)
    print(Y_pred)

    if(prefix=="train_eval"):
        heading="Model evaluation for Train data"
        color="red"
    elif(prefix=="test_eval"):
        heading="Model evaluation for Test data"
        color="green"
    elif(prefix=="validation" or re.search("validation",prefix)):
        heading="Model evaluation for Validation data"
        color="black"

    
    outfileHTML.write("<h2 style=text-align:center;color:"+color+">"+"----------------"+heading+"--------------------"+"</h2>"+"\n")
    outfileHTML.write("<h3>"+"Min Square Error for the Model"+"</h3>"+"\n")
    outfileHTML.write("<h4>"+"MSE of "+prefix+" set:{}".format(metrics.mean_squared_error(Y_test, Y_pred))+"</h4>")

    #Converting numpy.ndarray to dataframes
    column_dict_Y_test=dict()
    column_dict_Y_predict=dict()
    Y_test_df = pd.DataFrame(data=Y_test, columns=select_label_var_list)
    Y_pred_df = pd.DataFrame(data=Y_pred, columns=select_label_var_list)


    ##CALCULATE the RATIO of RMSE to Orignal Stdev for each label-feature
    ratio_low_dict=dict()
    ratio_high_dict=dict()
    ratio_dict=dict()
    mean_dict=dict()
    rmse_dict=dict()
    r2_score_dict=dict()
    std_dict=dict()
    n_rows=(len(select_label_var_list)/5)+1
    n_cols=5
    n_plots=len(select_label_var_list)
    print("performing further calculations for "+prefix)
    for i in select_label_var_list:
        Y_test_c=pd.DataFrame.to_numpy(Y_test_df[[i]])
        Y_pred_c=pd.DataFrame.to_numpy(Y_pred_df[[i]])
        rmse=sqrt(metrics.mean_squared_error(Y_test_c, Y_pred_c))
        r2score=r2_score(Y_test_c, Y_pred_c)
        mean=np.mean(Y_test_c)
        stdev=np.std(Y_test_c)
        if(stdev==0):
            if(mean==0):
                print(" The feature column "+i+" has mean and stdev values as zero. Its rmse is "+ str(rmse) + ". Not considering it for rmse:stdev ratio calculation")
                ratio="NA"
            else:
                if(rmse!=0):
                    print("Note:The feature column "+i+" has a non-zero mean and stdev is zero. It seems the feature column is monotonic. Its rmse is "+str(rmse)+" .Not considering the feature for rmse:stdev ratio calculation")
                    ratio="NA"
                elif(rmse==0):
                    print("Note:The feature column "+i+" has mean, stdev and rmse as zero. Hence excluding it from rmse:stdev ratio calculation")
                    ratio="NA"
        else:
            ratio=abs(rmse)/abs(stdev)
        mean_dict.update({i:mean})
        rmse_dict.update({i:rmse})
        r2_score_dict.update({i:r2score})
        std_dict.update({i:stdev})
        if(ratio=="NA"):
            ratio=-1.0
        elif(ratio <= 1.0):
            ratio_low_dict.update({i:ratio})
        elif(ratio > 1.0):
            ratio_high_dict.update({i:ratio})
        ratio_dict.update({i:ratio})
    label_header_low_ratio=list(ratio_low_dict.keys())


    mean_header=list(mean_dict.keys())
    mean_df=pd.DataFrame.from_dict(mean_dict,orient='index',columns=['Observed Mean'])
    
    std_header=list(std_dict.keys())
    std_df=pd.DataFrame.from_dict(std_dict,orient='index',columns=['Observed Stdev'])

    rmse_header=list(rmse_dict.keys())
    rmse_df=pd.DataFrame.from_dict(rmse_dict,orient='index',columns=['RMSE between observed and predicted values'])

    ratio_header=list(ratio_dict.keys())
    ratio_df=pd.DataFrame.from_dict(ratio_dict,orient='index',columns=['Ratio_of_RMSE_and_Stdev'])

    r2_score_header=list(r2_score_dict.keys())
    r2_score_df=pd.DataFrame.from_dict(r2_score_dict,orient='index',columns=['r2_score'])

    rmse_n_mean_df = rmse_df.merge(mean_df, how='outer', left_index=True, right_index=True)
    rmse_n_mean_n_std_df = rmse_n_mean_df.merge(std_df,how='outer', left_index=True, right_index=True)
    rmse_n_mean_n_std_n_ratio_df = rmse_n_mean_n_std_df.merge(ratio_df,how='outer', left_index=True, right_index=True)
    rmse_n_mean_n_std_n_ratio_n_r2_score_df = rmse_n_mean_n_std_n_ratio_df.merge(r2_score_df,how='outer', left_index=True, right_index=True)
    
    count=0
    Only_ratio_n_mean_df=rmse_n_mean_n_std_n_ratio_df[["Observed Mean","Ratio_of_RMSE_and_Stdev"]]
    Only_ratio_n_mean_df.plot.bar()
    rmse_n_mean_n_std_n_ratio_df.to_csv("/data/"+tagDir+prefix+'_'+model_type+'_rmse_mean_std_and_ratio.csv')
    rmse_n_mean_n_std_n_ratio_n_r2_score_df.to_csv("/data/"+tagDir+prefix+'_'+model_type+'_rmse_mean_std_and_ratio_and_r2_score.csv')
    mp.xticks(rotation=90)
    mp.ylabel('Ratio_of_RMSE_and_Stdev')
    mp.xlabel(label_type+" Features")
    mp.title("Observed Mean and Ratio_of_RMSE_and_Stdev",size=12)
    mp.legend(loc=(1.04,0.5))
    mp.savefig("/data/"+tagDir+prefix+'_'+model_type+'_rmse_mean_std_and_ratio.png',bbox_inches='tight')
    mp.clf()

    Only_ratio_n_r2_score_df=rmse_n_mean_n_std_n_ratio_n_r2_score_df[["Ratio_of_RMSE_and_Stdev","r2_score"]]
    Only_ratio_n_r2_score_df.plot.bar()
    mp.xticks(rotation=90)
    mp.ylabel('r2_score')
    mp.xlabel(label_type+" Features")
    mp.title("Ratio_of_RMSE_and_Stdev and r2_score",size=12)
    mp.legend(loc=(1.04,0.5))
    mp.savefig("/data/"+tagDir+prefix+'_'+model_type+'_rmse_std_and_ratio_and_r2_score.png',bbox_inches='tight')
    mp.clf()
    Y_test_low_ratio_df=Y_test_df[label_header_low_ratio]
    Y_pred_low_ratio_df=Y_pred_df[label_header_low_ratio]
    count=0
    for c in Y_pred_low_ratio_df.columns:
        mp.scatter(Y_test_low_ratio_df[c], Y_pred_low_ratio_df[c], label=c, marker=count)
        count=count+1
        if(count==11):
            count=0
    mp.xlabel('Actual_values')
    mp.ylabel('Predicted_values')
    mp.title("Actual_values v/s Predicted Values - for features with Low RMSE:Actual_Stdev",size=9)
    if len(label_header_low_ratio) <= 40:
        mp.legend(loc=(1.04,0))
    else:
        mp.legend(bbox_to_anchor=(1.04, 1.04, 2.04, 2.04), loc='upper left', ncol=2, mode="expand")
    mp.savefig("/data/"+tagDir+prefix+'_'+model_type+'_test_result_plots_low_ratio.png',bbox_inches='tight')
    mp.clf()

    label_header_high_ratio=list(ratio_high_dict.keys())
    
    Y_test_high_ratio_df=Y_test_df[label_header_high_ratio]
    Y_pred_high_ratio_df=Y_pred_df[label_header_high_ratio]
    
    count=0

    for d in Y_pred_high_ratio_df.columns:
        mp.scatter(Y_test_high_ratio_df[d], Y_pred_high_ratio_df[d], label=d, marker=count)
        count=count+1
        if(count==11):
            count=0
    mp.xlabel('Actual_values')
    mp.ylabel('Predicted_values')
    mp.title("Actual_values v/s Predicted Values - for features with high RMSE:Actual_Stdev",size=9)
    mp.legend(loc=(1.04,0))
    mp.savefig("/data/"+tagDir+prefix+'_'+model_type+'_test_result_plots_high_ratio.png',bbox_inches='tight')
    mp.clf()

    

    #Merging mean_df and rmse_dict
    if(len(ratio_low_dict)>1):
        ratio_low_df=pd.DataFrame.from_dict(ratio_low_dict,orient='index',columns=['RMSE/Stdev'])
        ratio_low_df.to_csv("/data/"+tagDir+prefix+"_"+model_type+"_Labels_with_Low_Ratio.csv")
        ratio_low_df.plot(kind='bar')
        mp.ylabel('RMSE/Stdev')
        mp.xlabel('Label Features')
        mp.xticks(rotation=90)
        if len(label_header_low_ratio) <= 40:
            mp.xticks(size=5)
        else:
            mp.xticks(size=3)
        mp.title("Low RMSE/Stdev for the label features",size=12)
        mp.savefig("/data/"+tagDir+prefix+'_'+model_type+'_Low_Ratio_plot.png',orientation='landscape',dpi=100,bbox_inches='tight')
        mp.clf()
    else:
        print("Only 1 key:value pair in ratio_low_dict, so not proceeding with its plotting")
    if(len(ratio_high_dict)>1):
        ratio_high_df=pd.DataFrame.from_dict(ratio_high_dict,orient='index',columns=['RMSE/Stdev'])
        ratio_high_df.to_csv("/data/"+tagDir+prefix+"_"+model_type+"_Labels_with_High_Ratio.csv")
        mp.plot(ratio_high_df)
        mp.ylabel('RMSE/Stdev')
        mp.xlabel('Label Features')
        mp.xticks(rotation=90)
        mp.xticks(size=4)
        mp.title("High RMSE/Stdev for the label features",size=12)
        mp.savefig("/data/"+tagDir+prefix+'_'+model_type+'_High_Ratio_plot.png',orientation='landscape',dpi=100,bbox_inches='tight')
        mp.clf()
    else:
        print("Only 1 key:value pair in ratio_high_dict, so not proceeding with its plotting")

    if(prefix=="train_eval"):
        heading="Model evaluation for Train data"
    elif(prefix=="test_eval"):
        heading="Model evaluation for Test data"
    elif(prefix=="validation"):
        heading="Model evaluation for Validation data"

    outfileHTML.write("<h4 stype=text-align:center;color:brown>"+"No. of features showing LOW 'RMSE/Stdev' (<=1.0): "+"\n"+str(len(label_header_low_ratio))+"</h4>")
    outfileHTML.write("<h5>"+"All such features with their Low 'RMSE/Stdev' values could be found in output file: "+prefix+"_"+model_type+"_Labels_with_Low_Ratio.csv"+"</h4>"+"\n")
    outfileHTML.write("<h4 stype=text-align:center;color:brown>"+"No. of features showing HIGH 'RMSE/Stdev' (>1.0): "+"\n"+str(len(label_header_high_ratio))+"</h4>")
    outfileHTML.write("<h5>"+"All such features with their High 'RMSE/Stdev' values could be found in output file: "+prefix+"_"+model_type+"_Labels_with_High_Ratio.csv"+"</h4>"+"\n"+"\n")
    outfileHTML.write("<h3>"+heading+" for label features showing Low 'RMSE/Stdev' (<=1.0)"+"</h3>"+"\n")
    
    data_image1 = open("/data/"+tagDir+prefix+'_'+model_type+'_test_result_plots_low_ratio.png', 'rb').read().encode('base64').replace('\n', '')
    img_tag1 = '<img src="data:image/png;base64,{0}">'.format(data_image1)
    outfileHTML.write(img_tag1+"\n")

    if(os.path.exists("/data/"+tagDir+prefix+'_'+model_type+'_Low_Ratio_plot.png')):
        data_image2 = open("/data/"+tagDir+prefix+'_'+model_type+'_Low_Ratio_plot.png', 'rb').read().encode('base64').replace('\n', '')
        img_tag2 = '<img src="data:image/png;base64,{0}">'.format(data_image2)
        outfileHTML.write(img_tag2+"\n")
    
    data_image3 = open("/data/"+tagDir+prefix+'_'+model_type+'_rmse_mean_std_and_ratio.png', 'rb').read().encode('base64').replace('\n', '')
    img_tag3 = '<img src="data:image/png;base64,{0}">'.format(data_image3)
    outfileHTML.write(img_tag3+"\n")

    
    outfileHTML.write('Content-type: text/html\n\n'+""+"\n"+'<link href="default.css" rel="stylesheet" type="text/css" />')

    decision_thresholds=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

    ##Adding feature-importance calculations
    FI_handle=open("/data/"+tagDir+prefix+"_"+"FeatureImportances.txt",'w+')
    X_features=list(test.columns)
    Y_features=list(Y_test.columns)
    if(gridsearch == 'True'):
        FI_handle.close()
    elif(model_type!="DecisionTree"):
        importance = model.coef_
        FI_handle.write("Y_Features, X_features_scored:"+str(X_features)+"\n")
        for u,k in enumerate(importance):
            FI_handle.write(str(Y_features[u])+","+str(k)+"\n")
        FI_handle.close()
    else:
        importance = model.feature_importances_
        FI_handle.write("Features, Score"+"\n")
        for u,k in enumerate(importance):
            FI_handle.write(str(X_features[u])+","+str(k)+"\n")
        FI_handle.close()
    


    ##CALCULATE AUCs BASED on the decision thresholds list above.
    AUC_fH=open("/data/"+tagDir+prefix+"_"+"AUC_values.txt",'w+')
    AUC_fH.write("label"+"\t"+"AUC_value"+"\t"+"decision_threshold"+"\n")
    label_resulted=[]
    for l in select_label_var_list:
        print l
        decision_threshold_list=[]
        AUC_value_list=[]
        for i in decision_thresholds:
            print i
            Y_test_df_n=binarize(pd.DataFrame.to_numpy(abs(Y_test_df[[l]])),threshold=i)
            Y_pred_df_n=Y_pred_df[[l]]##Keeping Ypredict to be continuous,i.e. as it is.
            if(np.all((Y_test_df_n==0))):
                print("For decision threshold "+str(i)+":")
                print("Seems all values for label column "+l+" are zero. Hence not considering it for decision_threshold vs AUC plot")
            else:
                fpr, tpr, thresholds = roc_curve(Y_test_df_n,Y_pred_df_n)
                print thresholds
                print fpr
                print tpr
                AUC_value=auc(fpr,tpr)
                if(AUC_value>0.9 and r2_score_dict[l]>0.25):
                    if l not in label_resulted:
                        label_resulted.append(l)
                AUC_fH.write(l+"\t"+str(AUC_value)+"\t"+str(i)+"\n")
                if(math.isnan(float(AUC_value))):
                    print "###AUC_value is nan#######"
                    print fpr; print tpr; print "#################"
                    continue
                AUC_value_list.append(AUC_value)
                decision_threshold_list.append(i)
        if(len(AUC_value_list)>=8 and len(decision_threshold_list)>=8):
            mp.plot(decision_threshold_list, AUC_value_list, label = '%s' % (l), linewidth=1, alpha=3)
            #label_resulted.append(l)
        else:
            print l+" has value list is not greater than or equal to 8"
    AUC_fH.close()

    mp.legend(loc=(1.04, 0))
    mp.title(prefix+" AUC_for_Decision_thresholds",size=12)
    mp.xlabel('absolute decision thresholds')
    mp.ylabel('AUC')
    mp.ylim(bottom=0.0)
    mp.xlim(left=0.1, right=0.9)
    mp.savefig("/data/"+tagDir+prefix+'_'+model_type+'_AUC_for_decision_thresholds.png',orientation='landscape',dpi=100,bbox_inches='tight')
    mp.clf()
    data_image4 = open("/data/"+tagDir+prefix+'_'+model_type+'_AUC_for_decision_thresholds.png', 'rb').read().encode('base64').replace('\n', '')
    img_tag4 = '<img src="data:image/png;base64,{0}">'.format(data_image4)
    outfileHTML.write(img_tag4+"\n")

    outfileHTML.write('Content-type: text/html\n\n'+""+"\n"+'<link href="default.css" rel="stylesheet" type="text/css" />')
    return label_resulted

## MODEL PREDICTION function    
def predict(model , test):
    '''
        Get model predictions
        Returns predicted data as a numpy array

        :param model: Model used to perform prediction
        :type model: model
        :param test: data to perform prediction on
        :type test: dataframe
    '''
    Y_pred = model.predict(test)
    return Y_pred

##PROCESS FUNCTION       
def process(data_, label_, data_type, label_type, corr_method, corr_threshold, pVal_adjust_method, data_normalize_method, label_normalize_method, cv_par, scoring_par, mode, model_type, load_model, params, grid_search, param_grid, prediction_out, select_label_headers_for_predict, select_data_headers_for_predict, featureSelFrmModel_flag=None):
    '''
        PROCESS FUNCTION: Entailing the entire process starting from training through testing and validation per the modes of operation specified by user.
        Returns None

        :param data_: data file.
        :type data_: file
        :param label_: label file.
        :type label_: file
        :param data_type: A string indicating the type of data, ex: "imaging" or "gene".
        :type data_type: str
        :param label_type: A string indicating the type of label, ex: "imaging" or "gene".
        :type label_type: str
        :param corr_method: Correlation method used, options: {'Pearson', 'Spearman'}
        :type corr_method: str
        :param corr_threshold: Threshold for correlation co-efficient, options: {-1.0,0.0,0.1...1.0}, default:0.5
                            ##Note: To silent correlation threshold filtering, set the threshold to -1.0
        :type corr_threshold: float
        :param pVal_adjust_method: The method to correct p-values, options: {'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'}, default: BH
        :type pVal_adjust_method: str
        :param data_normalize_method: Normalization method for data, options: {'min_max' , 'Stand_scaler' , 'MaxAbsScaler', 'none'}, default: None
        :type data_normalize_method: str
        :param label_normalize_method: Normalization method for label, options: {'min_max' , 'Stand_scaler' , 'MaxAbsScaler', 'none'}, default: None
        :type label_normalize_method: str
        :param cv_par: K-fold Cross validation splitter number, default: 2
        :type cv_par: int
        :param scoring_par: Score metric, default: neg_mean_square_error
        :type scoring_par: str
        :param mode: Mode of run, options: {'Train','predict','validate'}. ##Note: This parameter is no more in this function and would be removed from argument list in next version.
        :type mode: str
        :param model_type: Model type, options: {'DecisionTree','LinearRegression', 'LinearModel' , 'LASSO', 'multiTaskLASSO', 'multiTaskLinearModel'}, default: None
        :type model_type: str
        :param load_model: Model pickel (.pkl) file to load in case mode set to "predict" or "validate"
        :type load_model: file
        :param params: Model parameters. If empty, default parameters used per scikit-learn library.
        :type params: dict
        :param gridsearch: Indicating whether to perform gridsearch for training the model, options: 'True' or 'False', default: 'False'
        :type gridsearch: str
        :param param_grid: Parameters for gridsearch for the given model type.
        :type param_grid: dict
        :param prediction_out: Output file to write predicted output to, in case mode set to "predict".
        :type prediction_out: file
        :param select_data_headers_for_predict: List of data features to used when mode is set to "predict" or "validate".
        :type select_data_headers_for_predict: list
        :param select_label_headers_for_predict: List of data features to used when mode is set to "predict" or "validate".
        :type select_label_headers_for_predict: list
        :param featureSelFrmModel_flag: flag to run 'FeatureSelection' mode.
        :type featureSelFrmModel_flag: bool
    '''
    if(featureSelFrmModel_flag==None or featureSelFrmModel_flag==0):
        tagDir=""
    elif(featureSelFrmModel_flag==1):
        tagDir="FeaturesSelFrmModel.txt/"
    if os.path.isfile(data_):
        if os.path.getsize(data_)!=0:
            dataframe = read_dataset(data_)
        else:
            sys.exit("Size of data_ file"+data_+" is zero")
    else:
        sys.exit("data_ file:"+data_+" does not exists as a regular file")
    if mode != 'predict':
        if(os.path.isfile(label_)):
            if(os.path.getsize(label_)!=0):
                label = read_dataset(label_)
            else:
                sys.exit("Size of label_ file:"+label_+"  is zero")
        else:
            sys.exit("data_ file:"+data_+" does not exists as a regular file")
    else:
        label='NA'
    dataframe , label, sampleIDs, label_header, dataframe_header =  preprocessing(dataframe , label, data_type, label_type, mode, tagDir)
    ##When mode selected is "Train"
    if(mode!='predict' and mode!='validate'):
        if(featureSelFrmModel_flag==0 or featureSelFrmModel_flag==None):
            select_data_var_list,select_label_var_list = correlation(dataframe,label,dataframe_header,label_header, model_type, corr_method, corr_threshold, pVal_adjust_method, tagDir)
        elif(featureSelFrmModel_flag==1):
            select_data_var_list=dataframe_header; select_label_var_list=label_header
        print select_label_var_list
        outfileHTML=open("/data/"+tagDir+model_type+".output.html",'a')
        outfileHTML.write("<h1 style=text-align:center>"+"----------------------------------Features with highly significant correlations-------------------------------------"+"</h1>"+"\n")
        dataframe = dataframe[select_data_var_list]
        outfileHTML.write("<h2>"+"Below is the list of "+data_type+" features"+"</h2>"+"\n")
        outfileHTML.write("<h5>"+str(list(dataframe.keys()))+"</h5>"+"\n")
        label = label[select_label_var_list]
        print label
        print dataframe

        outfileHTML.write("<h2>"+"Below is the list of "+label_type+" features"+"</h2>"+"\n")
        outfileHTML.write("<h5>"+str(list(label.keys()))+"</h5>"+"\n")
        outfileHTML.write("<h2>"+"Number of "+data_type+" features"+"</h2>"+"\n")
        outfileHTML.write("<h3>"+str(len(select_data_var_list))+"</h3>"+"\n")
        outfileHTML.write("<h2>"+"Number of "+label_type+" features"+"</h2>"+"\n")
        outfileHTML.write("<h3>"+str(len(select_label_var_list))+"</h3>"+"\n")
        outfileHTML.close()
    ##When mode selected is "predict"
    elif(mode=='predict'):
        if(len(select_data_headers_for_predict)!=0):
            dataframe = dataframe[select_data_headers_for_predict]
        if(len(select_label_headers_for_predict)!=0):
            label_header_for_predict = select_label_headers_for_predict
    ##When mode set to "validate"
    elif(mode=='validate'):
        if(len(select_data_headers_for_predict)!=0):
            dataframe = dataframe[select_data_headers_for_predict]
        if(len(select_label_headers_for_predict)!=0):
            select_label_var_list_for_validate=select_label_headers_for_predict
            label = label[select_label_var_list_for_validate]
        else:
            select_label_var_list_for_validate=label_header

    ##TRAINING and TESTING the model
    if mode == 'Train':
        train , Y_train ,test , Y_test = splitdata(dataframe , label, test_size, mode, data_normalize_method, label_normalize_method, data_type, label_type, select_data_var_list, select_label_var_list, tagDir)
        print("Staring Training of :{}".format(model_type))
        model = BuildModel(train , Y_train , test , Y_test , model_type, params, cv_par, scoring_par, grid_search, param_grid, select_label_var_list, select_data_var_list, data_type, label_type, featureSelFrmModel_flag, tagDir, trainmodel = 'True')
        if save == 'True':
            try:
                joblib.dump(model, str(save_dir) + str(tagDir) + str(model_type)+ ".pkl" )
                print("Model saved at:{}".format(str(save_dir) + str(tagDir) + str(model_type)+ ".pkl"))
            except OSError:
                print("Saving model failed")
        else:
            print("Please provide save dir")
    ##PREDICTION
    elif mode == 'predict':
        print("Performing Prediction")
        try:
            model = joblib.load(str(load_model))
            print(model.get_params)

            outfileHTML=open("/data/"+tagDir+model_type+".output.html",'a')
            store_params=model.get_params();
            outfileHTML.write("<h2 style=text-align:center;color:blue>"+"------------------------Model Summary-----------------------"+"</h2>")
            outfileHTML.write("<h3>"+"Model Parameters:"+"</h3>"+"\n")
            for i in store_params.keys():
                outfileHTML.write("<h4>"+str(i)+":"+str(store_params[i])+"</h4>")
            outfileHTML.write("<h2 style=text-align:center;color:green>"+"------------------------Samples for Prediction-----------------------"+"</h2>")
            outfileHTML.write("<h3>"+"No. of samples input for prediction: "+str(len(sampleIDs))+"</h3>")

            ##Normalizing data
            outfileHTML.write("<h3>"+"performing "+data_normalize_method+" normalization for "+data_type+" features for prediction set"+"</h3>"+"\n")
            if(len(select_data_headers_for_predict)!=0):
                dataframe = normal_dataframe(dataframe, data_normalize_method, select_data_headers_for_predict)
            else:
                dataframe = normal_dataframe(dataframe, data_normalize_method, dataframe_header)

            Y_pred = predict(model, dataframe)
            if prediction_out == "NULL":
               prediction_out = "prediction_out.txt"
            outfileH=open("/data/"+tagDir+prediction_out,'w');

            outfileHTML.write("<h3>"+"All predicted values for labels available in the file: "+prediction_out+"</h3>")

            if(label_header_for_predict):
                outfileH.write("sampleID"+"\t"+"\t".join(label_header_for_predict)+"\n")
                outfileH.close()
                outfileH=open("/data/"+tagDir+prediction_out,'a')
            if len(sampleIDs) != len(Y_pred):
                outfileHTML.write("<h5 style=color:red>"+"No. of samples provided as input DO NOT MATCH WITH No. of samples predicted. Therefore, no prediction performed. Kindly investigate the log file for errors"+"</h5>")
                sys.exit("The number of samples in ID column in the data does not match with the number of samples for which predicted values were obtained. Kindly check your data file for possible issues.")
            else:
                outfileHTML.write("<h5 style=color:green>"+"No. of samples provided as input MATCH WITH No. of samples predicted"+"</h5>")
            outfileHTML.close()
            count=0
            for i in Y_pred:
                outfileH.write(sampleIDs[count]+"\t")
                k=len(i)
                for j in range(0,k):
                    if j != k-1:
                        outfileH.write(str(i[j])+"\t")
                    else:
                        outfileH.write(str(i[j]))
                outfileH.write("\n")
                count=count+1
            outfileH.close()
        except IOError:
            print("Probably saved model file is not provided/check path")
        except (IndexError, ValueError), err_:
            print("Seems an IndexError or a ValueError was encountered during prediction. Error msg is as follows: "+str(err_))
    ##VALIDATION of the model
    elif mode == 'validate':
        print("Performing validation")
        try:
            model = joblib.load(str(load_model))
            print(model.get_params)
            outfileHTML=open("/data/"+tagDir+model_type+".output.html",'a')
            store_params=model.get_params();
            outfileHTML.write("<h2 style=text-align:center;color:blue>"+"------------------------Model Summary-----------------------"+"</h2>")
            outfileHTML.write("<h3>"+"Model Parameters:"+"</h3>"+"\n")
            for i in store_params.keys():
                outfileHTML.write("<h4>"+str(i)+":"+str(store_params[i])+"</h4>")
            outfileHTML.write("<h2 style=text-align:center;color:green>"+"------------------------Samples for Validation-----------------------"+"</h2>")
            outfileHTML.write("<h3>"+"No. of samples used for validation: "+str(len(sampleIDs))+"</h3>")
            
            ##Normalizing data
            outfileHTML.write("<h3>"+"performing "+data_normalize_method+" normalization for "+data_type+" features"+"</h3>"+"\n")
            if(len(select_data_headers_for_predict)!=0):
                dataframe = normal_dataframe(dataframe, data_normalize_method, select_data_headers_for_predict)
            else:
                dataframe = normal_dataframe(dataframe, data_normalize_method, dataframe_header)

            ##Normalizing label
            outfileHTML.write("<h3>"+"performing "+label_normalize_method+" normalization for "+label_type+" features"+"</h3>"+"\n")
            if(len(select_label_var_list_for_validate)!=0):
                label = normal_dataframe(label, label_normalize_method, select_label_var_list_for_validate)
            else:
                label = normal_dataframe(label, label_normalize_method, label_header)
            outfileHTML.close()
            
            ##Evaluate
            evaluate(model, dataframe , label, select_label_var_list_for_validate, 'validation', model_type, data_type, label_type, tagDir)
        except IOError:
            print("Probably saved model file is not provided/check path")
    else:
        sys.exit("Invalid mode. Should be either Train, predict or validate. Please check your config file.")

## MAIN function. The main starting point of the script. Scans parameters from config.ini, defines variable and passes on to the process function called below.
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data' , help = "data file" , default = 'NULL', type = str)
    parser.add_argument('--label', help = "label file", default = 'NULL', type = str)
    parser.add_argument('--config', help = "config file", default = 'NULL', type = str)
    parser.add_argument('--model_file', help = "model file, .pkl type", default = 'NULL', type = str)
    parser.add_argument('--prediction_out', help = "prediction output", default = 'NULL', type = str)
    args = parser.parse_args()
    
    #Parsing config
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(args.config)
    save_dir = config['paths']['save_dir']
    model_type = config['modes']['model']
    mode = config['modes']['mode']
    if(mode=="predict" or mode=="validate"):
        load_model=args.model_file
        model_name_extract=load_model.split('.')[0].split('/')[-1]
        if(model_type!=model_name_extract):
            sys.exit("The model_type specified in config file does not match with model_file loaded")
    elif(mode=="Train"):
        load_model="NULL"
    else:
        sys.exit("Invalid mode. Should be either Train, predict or validate. Please check your config file.")
    data_type = config['modes']['data_type']
    label_type = config['modes']['label_type']
    data_normalize_method = config['modes']['data_normalize_method']
    label_normalize_method = config['modes']['label_normalize_method']
    save = config['modes']['save']
    grid_search = config['modes']['grid_search']
    cv_par = int(config['modes']['cv'])
    scoring_par = config['modes']['scoring']
    test_size = float(config['modes']['test_size'])
    corr_method = config['modes']['correlation_method']
    corr_threshold = float(config['modes']['correlation_threshold'])
    pVal_adjust_method = config['modes']['pVal_adjust_method']

    ##getting data and label features from respective parameters in 'modes' section, useful when validating a model that ended up trained on selected features.
    select_data_headers_for_predict = config['modes']['select_data_headers_for_predict']
    select_label_headers_for_predict = config['modes']['select_label_headers_for_predict']
    if isinstance(select_data_headers_for_predict,unicode):
        try:
            select_data_headers_for_predict=ast.literal_eval(select_data_headers_for_predict)
        except ValueError:
            select_data_headers_for_predict=select_data_headers_for_predict.encode('utf-8')
    if isinstance(select_label_headers_for_predict,unicode):
        try:
            select_label_headers_for_predict=ast.literal_eval(select_label_headers_for_predict)
        except ValueError:
            select_label_headers_for_predict=select_label_headers_for_predict.encode('utf-8')
    
    #storing params for model from sections in config
    section_flag=0; params = {}; param_grid = {};
    for each_section in config.sections():
        #print each_section
        if each_section == model_type and grid_search == 'False':
            print "Section for model "+each_section+" detected."
        elif each_section == model_type+"_grid" and grid_search == 'True':
            section_flag=1
            print "Grid section "+each_section+" detected and grid_search is set True."
        else:
            continue
        for (each_key, each_val) in config.items(each_section):
            #print each_key; print each_val
            if each_val.replace('.','',1).isdigit() == True:
                if each_val.isdigit():
                    each_val = int(each_val)
                else:
                    each_val = float(each_val) 
            elif isinstance(each_val,str):
                if each_val == 'True':
                    each_val = True
                elif each_val == 'False':
                    each_val = False
                elif each_val == 'None':
                    each_val = None
            elif isinstance(each_val,unicode):
                try:
                    each_val=ast.literal_eval(each_val)
                except ValueError:
                    each_val=each_val.encode('utf-8')
            if each_val != '':
                if section_flag == 0 :
                    params[each_key] = each_val
                    #print str(params)
                elif section_flag == 1 :
                    param_grid[each_key] = each_val

    if len(params) == 0:
        print_params = "default"
    else:
        print_params = params
    if len(param_grid) == 0:
        print_param_grid = "default"
    else:
        print_param_grid = param_grid
    outfileHTML=open("/data/"+model_type+".output.html",'w')
    FeaturesSelFrmModel_dir_path="/data/"+"FeaturesSelFrmModel.txt/"
    if not os.path.exists(FeaturesSelFrmModel_dir_path):
        os.makedirs(FeaturesSelFrmModel_dir_path)## Making FeaturesSelFrmModel dir
    outfileHTML2=open("/data/"+"FeaturesSelFrmModel.txt/"+model_type+".output.html",'w')##Getting it ready for FeatureSelFrmModel mode as well.
    outfileHTML.write("<h1 style=text-align:center;color:red;>"+"Radiogenomics Analysis Report"+"</h1>"+"\n")
    outfileHTML2.write("<h1 style=text-align:center;color:red;>"+"Radiogenomics Analysis Report with Features_Selection_From_Model mode activated"+"</h1>"+"\n")
    #write date and time of the report
    from datetime import datetime
    datetime_now = datetime.now()
    datetime_string = datetime_now.strftime("%d/%m/%Y %H:%M:%S")
    outfileHTML.write("<h5 style=text-align:center;>"+datetime_string+"</h5>"+"\n")
    outfileHTML.write("<h2 style=text-align:center;color:brown>"+"----------------------------Model inputs-------------------------------"+"</h2>")
    outfileHTML.write("<h3>"+"Mode:{}".format(mode)+"</h3>")
    outfileHTML.write("<h3>"+"Model:{}".format(model_type)+"</h3>")
    outfileHTML.write("<h3>"+"Params:{}".format(print_params)+"</h3>")
    outfileHTML2.write("<h5 style=text-align:center;>"+datetime_string+"</h5>"+"\n")
    outfileHTML2.write("<h2 style=text-align:center;color:brown>"+"----------------------------Model inputs-------------------------------"+"</h2>")
    outfileHTML2.write("<h3>"+"Mode:{}".format(mode)+"</h3>")
    outfileHTML2.write("<h3>"+"Model:{}".format(model_type)+"</h3>")
    outfileHTML2.write("<h3>"+"Params:{}".format(print_params)+"</h3>")
    if grid_search == "True":
        outfileHTML.write("<h3>"+"Grid_Params:{}".format(print_param_grid)+"</h3>")
        outfileHTML2.write("<h3>"+"Grid_Params:{}".format(print_param_grid)+"</h3>")
    outfileHTML.close()
    outfileHTML2.close()
    
    ##CALL the PROCESS FUNCTION.
    process(args.data, args.label, data_type, label_type, corr_method, corr_threshold, pVal_adjust_method, data_normalize_method, label_normalize_method, cv_par, scoring_par, mode, model_type, load_model, params, grid_search, param_grid, args.prediction_out, select_label_headers_for_predict, select_data_headers_for_predict, featureSelFrmModel_flag=0)
    ##Calling process again with feature selection flag (featureSelFrmModel_flag) equal to 1, so that the feature selection performed using SelectFromModel module (i.e., Feature Selection using Feature Importances from model), and the correlation module is skipped completely.
    if mode == "Train" and grid_search == "False":
        process(args.data, args.label, data_type, label_type, corr_method, corr_threshold, pVal_adjust_method, data_normalize_method, label_normalize_method, cv_par, scoring_par, mode, model_type, load_model, params, grid_search, param_grid, args.prediction_out, select_label_headers_for_predict, select_data_headers_for_predict, featureSelFrmModel_flag=1)
