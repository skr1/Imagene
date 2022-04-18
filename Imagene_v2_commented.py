##Imagene.py
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
##Author: Shrey Sukhadia
#!/usr/bin/python
from cProfile import label
import matplotlib as mpl
mpl.use('Agg')
import os, re, sys, math
import argparse
import numpy as np
import pandas as pd
#from pandas.stats.api import ols
import joblib
import pickle
import configparser
import ast
import argparse
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
#model packages
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
#import statsmodels as statsm
#import statsmodels.api as smapi
#from mpl_toolkits.mplot3d import Axes3D
#import rpy2.robjects as R
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
##Correlation function establishes correlations between gene and imaging features and calculates FDR adjusted pvalues further to establish statistical significance of those correlations.
def correlation(d_1,d_2,d_1_header,d_2_header, model_type, corr_method, corr_threshold,pVal_adjust_method):
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
    #print(pcorr_PAJ_df_sorted.head(30))
    Master_df_sorted.to_csv("All_correlations.csv")

    Master_df_sorted_sgn=Master_df_sorted[Master_df_sorted['p_adjust']<0.05]
    Master_df_sorted_sgn_filtered=Master_df_sorted_sgn[abs(Master_df_sorted_sgn['correlation'])>corr_threshold]
    
    
    
    Master_df_sorted_sgn.to_csv("Significant_correlations.csv")
    Master_df_sorted_sgn_filtered.to_csv("Significant_correlations_gt_0.5.csv")
    #print(Master_df_sorted_sgn)
    List_of_Var1=Master_df_sorted_sgn_filtered["Var1"].tolist()
    List_of_Var2=Master_df_sorted_sgn_filtered["Var2"].tolist()
    fC_List_Var1=[]
    for i in List_of_Var1:
        fC=i.split("_")[0]
        fC_List_Var1.append(fC)
    fC_uniq = sorted(set(fC_List_Var1))

    #outfileHTML=open(model_type+".output.html",'w')
    image_tag_list=[]
    #image_count=0

    for i in fC_uniq:
        print(i)
        #fC_df=Master_df_sorted_sgn[Master_df_sorted_sgn['Var1'].str.match(i, case=True, flags=0)].head(50)
        #fC_df=Master_df_sorted_sgn[Master_df_sorted_sgn['Var1'].str.contains(i, regex=False, na=False)].head(50)
        fC_df=Master_df_sorted_sgn_filtered[Master_df_sorted_sgn_filtered['Var1'].str.contains(i, regex=False, na=False)]
        fC_df_pivot=fC_df.pivot_table(index="Var1",columns="Var2",values="correlation",fill_value=0)
        #fC_df_pivot = pd.pivot_table(fC_df, index = "Var1", values = ["Var2", "spearman_correlation"]).stack().reset_index(level = 1) 
        print("For "+i+" features:")
        print(fC_df_pivot)
        try:
            sb.clustermap(fC_df_pivot)
        except ValueError as err:
            mp.clf()
            sb.heatmap(fC_df_pivot)
        
        mp.title("Top significant correlations (FDR_adjusted_pValue<0.05) for "+i+" features", size=16)
        mp.xticks(rotation=90)
        mp.yticks(size=7)
        mp.savefig('/data/'+'Correlation_for_'+i+'_features.png',orientation='landscape',dpi=90,bbox_inches='tight')
        data_image = open("/data/"+'Correlation_for_'+i+'_features.png', 'rb').read().encode('base64').replace('\n', '')
        image_tag_list.append('<img src="data:image/png;base64,{0}" style="max-width:50%;">'.format(data_image))
        mp.clf()

    outfileHTML=open("/data/"+model_type+".output.html",'a')
    #outfileHTML.write("<h1 style=text-align:center;color:red;>"+"Radiogenomics Analysis Report"+"</h1>"+"\n")
    outfileHTML.write("<h2 style=text-align:center;>"+"--------------------------Multivariate Correlations ("+corr_method+" based)-----------------------"+"</h2>"+"\n".join(image_tag_list)+"\n")
    outfileHTML.close()





    return(sorted(set(List_of_Var1)),sorted(set(List_of_Var2)))

##Reading datasets
def read_dataset(dataset):
    if os.path.isfile(dataset):
        raw_dataframe = pd.read_csv(dataset, sep=',')
    else:
        print("{} should be csv/tsv file ".format(dataset))
        sys.exit()
    dataframe = raw_dataframe.copy()
    print("Shape of dataframe:{}".format(dataframe.shape))
    return dataframe
##Normalization function
def normal_dataframe(dataframe, norm_type, header, normalize = True):
    if normalize:
        if norm_type =='min_max':
            min_max_scaler = MinMaxScaler()
            dataframe_scaled = min_max_scaler.fit_transform(dataframe)
            #return dataframe_scaled
        
        elif norm_type == 'Stand_scaler':
            scaler = StandardScaler()
            dataframe_scaled = scaler.fit_transform(dataframe)
            #return dataframe_scaled
        
        elif norm_type == 'zscore':
            dataframe_scaled = stats.zscore(dataframe)
            #return dataframe_scaled
        
        elif norm_type == 'MaxAbsScaler':
            max_abs_scaler = MaxAbsScaler()
            dataframe_scaled = max_abs_scaler.fit_transform(dataframe)
            #return dataframe_scaled
        else:
            print("Invalid normalization type/method detected. Skipping normalization. If you wish to normalize this dataset, then correct the normalization_method in the config file and rerun. Proceeding with no normalization")
            return dataframe
        dataframe_scaled=pd.DataFrame(data=dataframe_scaled,columns=header)
        return dataframe_scaled
    else:
        return dataframe
##Preprocessing of datasets    
def preprocessing(dataframe , label, data_type, label_type, data_normalize_method, label_normalize_method , mode, checkNA = True):
    #try:
        outfileHTML=open("/data/"+model_type+".output.html",'a')
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
        
        
        if data_normalize_method != 'none':
            print("performing "+data_normalize_method+" normalization for "+data_type+" features")
            outfileHTML.write("<h3>"+"performing "+data_normalize_method+" normalization for "+data_type+" features"+"</h3>"+"\n")
            dataframe = normal_dataframe(dataframe ,data_normalize_method, dataframe_header)
            #print(dataframe)
        if isinstance(label,pd.DataFrame):
            if label_normalize_method != 'none':
                #print(label)
                print("performing "+label_normalize_method+" for "+label_type+" features")
                outfileHTML.write("<h3>"+"performing "+label_normalize_method+" for "+label_type+" features"+"</h3>"+"\n")
                label = normal_dataframe(label ,label_normalize_method, label_header)
        outfileHTML.close()
        #c_=pd.concat([dataframe1, dataframe2], axis=1)
        dataframe = dataframe.loc[:, (dataframe!=0).any(axis=0)]; nan_value = float("NaN"); dataframe.replace("", nan_value, inplace=True); dataframe=dataframe.dropna()
        dataframe_header=list(dataframe.keys())
        print(dataframe.shape[1])
        if isinstance(label,pd.DataFrame):
            label = label.loc[:, (label!=0).any(axis=0)]; label.replace("", nan_value, inplace=True); label=label.dropna()
            label_header=list(label.keys())
            print(label.shape[1])
        return dataframe , label, sampleIDs, label_header, dataframe_header

##Train-test split function    
def splitdata(dataframe , label, t_size, mode_):
    outfileHTML=open("/data/"+model_type+".output.html",'a')
    train, test , Y_train , Y_test = train_test_split(dataframe, label , test_size = t_size)
    #if mode_=="Train":
    outfileHTML.write("<h1 style=text-align:center;color:purple>"+"-------------------------------Number of Samples for Training and Testing---------------------------------"+"</h1>"+"\n")
    outfileHTML.write("<h3>"+"No. of samples for training:{}".format(len(train))+"</h3>"+"\n")
    outfileHTML.write("<h3>"+"No. of samples for test:{}".format(len(test))+"</h3>"+"\n")
    outfileHTML.close()
    print("Trainig data:{} , Testing data:{} ".format(len(train) ,len(test)))
    #if mode_=="validate":
    #    print("Custom test data:{} ".format(len(test)))
    return train , Y_train , test , Y_test

##Build models, multiple options of modeltypes accepted from user per the method list below
def BuildModel(train , Y_train , test , Y_test , method, params, cv_par, scoring_par, gridsearch, param_grid, select_label_var_list, data_type, label_type, trainmodel):
    '''
    initializing model and training the model
    '''
    outfileHTML=open("/data/"+model_type+".output.html",'a')
    outfileHTML.write("<h2 style=text-align:center;color:blue>"+"--------------------------Model Summary-----------------------"+"</h2>"+"\n")
    if method in ['DecisionTree','LinearRegression', 'LinearModel' , 'LASSO', 'multiTaskLASSO', 'multiTaskLinearModel']:
        if gridsearch == 'True':
            #param_grid_keys=list(param_grid.keys())
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
            #for i in param_grid_keys:
            #    if re.search('cv',i):
            #        cv_grid = param_grid[i]
            #        del param_grid[i]
            #    elif re.search('scoring',i):
            #        scoring_grid = param_grid[i]
            #        del param_grid[i]
            #    else:
            #        continue
        if method == 'LASSO':
            if gridsearch == 'True':
                try:
                    print(" Starting grid search for LASSO")
                    model = GridSearchCV(linear_model.Lasso(), param_grid=param_grid,cv=cv_grid,scoring=scoring_grid)
                    #print(model.get_params)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                #model = linear_model.Lasso(alpha=alpha_lasso)
                model = linear_model.Lasso(**params)
                #print(model.get_params)
                
        if method == 'DecisionTree':
            if gridsearch == 'True':
                try:
                    print("starting grid search for Decision Tress")
                    model = GridSearchCV(tree.DecisionTreeRegressor(), param_grid=param_grid, scoring=scoring_grid, cv=cv_grid)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                model = tree.DecisionTreeRegressor(**params)
                #print(model.get_params)
                
        if method == 'LinearRegression':
            if gridsearch == 'True':
                try:
                    print("starting grid search for Linear Regression")
                    model = GridSearchCV(LinearRegression(), param_grid=param_grid, scoring=scoring_grid, cv=cv_grid)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                model = LinearRegression(**params)
                #print(model.get_params)

        
        if method == 'LinearModel': 
            if gridsearch == 'True':
                #random_state=param_grid['random_state']
                #del new_param_grid['random_state']
                try:
                    print("Starting grid search for ElasticNet")
                    model = GridSearchCV(ElasticNet(), param_grid=param_grid , cv = cv_grid , scoring=scoring_grid)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                model = ElasticNet(**params)
                #print(model.get_params)

        if method == 'multiTaskLASSO':
            if gridsearch == 'True':
                try:
                    print(" Starting grid search for MultiTaskLASSO")
                    model = GridSearchCV(linear_model.MultiTaskLasso(), param_grid=param_grid,cv=cv_grid,scoring=scoring_grid)
                    #print(model.get_params)
                except:
                    print("Grid search status:{}".format(grid_search))
            else:
                #model = linear_model.Lasso(alpha=alpha_lasso)
                model = linear_model.MultiTaskLasso(**params)
                #print(model.get_params)

        if method == 'multiTaskLinearModel': 
            if gridsearch == 'True':
                #random_state=param_grid['random_state']
                #del new_param_grid['random_state']
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
        #try:
            '''
            performing training
            '''
            #outfileHTML=open(model_type+".output.html",'a')
            #outfileHTML.write("<h1>"+"--------------------------Model Summary-----------------------"+"</h1>"+"\n")
            outfileHTML.write("<h3>"+"Model Type : "+method+"</h3>"+"\n")
            if gridsearch == 'True':
                
                grid_result = model.fit(train ,Y_train)
                outfileHTML.write("<h4>"+"Grid Search Metrics"+"</h4>"+"\n")
                #outfileHTML.write("<h3>"+"Best Parameters: "+"</h3>"+"\n")
                #outfileHTML.write("<h3>"+ "scoring="+str(scoring_grid)+"</h3>"+"\n")
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
            
            #Y_pred_train = model.predict(train)
            #print("MSE of Train set:{}".format(metrics.mean_squared_error(Y_train, Y_pred_train)))
            #Y_pred = model.predict(test)
            #print("MSE of Test set:{}".format(metrics.mean_squared_error(Y_test, Y_pred)))

            outfileHTML.close()
            evaluate(model,train,Y_train,select_label_var_list,'train_eval',model_type, data_type, label_type)
            returned_label_resulted=evaluate(model,test,Y_test,select_label_var_list,'test_eval',model_type, data_type, label_type)
            ##Permut each label coulmn of the test dataset 100 times and validate the model
            for l in returned_label_resulted:
                for n in range(1,21):
                    Y_test[[l]]=np.random.permutation(Y_test[[l]])
                    evaluate(model, test, Y_test , select_label_var_list, 'validation'+'_permut_'+str(n)+'_'+str(l), model_type, data_type, label_type)
        #except:
            #print(" Issue with model training")
        
    return model
##Evaluate the models for test or validation.
def evaluate(model , test , Y_test, select_label_var_list, prefix, model_type, data_type, label_type):
    '''
    evalauting the model preformance 
    '''
    outfileHTML=open("/data/"+model_type+".output.html",'a')
    Y_pred = model.predict(test)

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
    #for i in label_header:
    #    column_dict_Y_test.update({i: Y_test[:, label_header.index(i)]})
    #    column_dict_Y_predict.update({i: Y_predict[:, label_header.index(i)]})
    #Y_test_df=pd.DataFrame(column_dict_Y_test)
    #Y_pred_df=pd.DataFrame(column_dict_Y_predict)
    Y_test_df = pd.DataFrame(data=Y_test, columns=select_label_var_list)
    #print Y_test_df
    Y_pred_df = pd.DataFrame(data=Y_pred, columns=select_label_var_list)
    #print Y_pred_df


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
        #if(i=='module1'):
        #    print Y_test_c
        #    print Y_pred_c
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
        #print i+"\t"+str(mean)
        #rmse_n_mean_df.update()
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
            #ratio_dict.update({i:ratio})
            #mean_for_rmse_high_dict.update({i:mean})
            #print("HIGH RMSE for "+i+" :{}".format(rmse))
            #mse_high_file.write("HIGH MSE for "+i+" :{}".format(mse))
            #mse_high_file.write("\n")
        ratio_dict.update({i:ratio})
    label_header_low_ratio=list(ratio_low_dict.keys())


    mean_header=list(mean_dict.keys())
    mean_df=pd.DataFrame.from_dict(mean_dict,orient='index',columns=['Observed Mean'])
    #mean_df.to_csv(prefix+"_"+model_type+"_mean.csv")
    
    std_header=list(std_dict.keys())
    std_df=pd.DataFrame.from_dict(std_dict,orient='index',columns=['Observed Stdev'])
    #std_df.to_csv(prefix+"_"+model_type+"_std.csv")

    rmse_header=list(rmse_dict.keys())
    rmse_df=pd.DataFrame.from_dict(rmse_dict,orient='index',columns=['RMSE between observed and predicted values'])

    ratio_header=list(ratio_dict.keys())
    ratio_df=pd.DataFrame.from_dict(ratio_dict,orient='index',columns=['Ratio_of_RMSE_and_Stdev'])

    r2_score_header=list(r2_score_dict.keys())
    r2_score_df=pd.DataFrame.from_dict(r2_score_dict,orient='index',columns=['r2_score'])
    #rmse_df.to_csv(prefix+"_"+model_type+"_rmse.csv")
    
    #print rmse_header
    #print mean_header

    rmse_n_mean_df = rmse_df.merge(mean_df, how='outer', left_index=True, right_index=True)
    rmse_n_mean_n_std_df = rmse_n_mean_df.merge(std_df,how='outer', left_index=True, right_index=True)
    rmse_n_mean_n_std_n_ratio_df = rmse_n_mean_n_std_df.merge(ratio_df,how='outer', left_index=True, right_index=True)
    rmse_n_mean_n_std_n_ratio_n_r2_score_df = rmse_n_mean_n_std_n_ratio_df.merge(r2_score_df,how='outer', left_index=True, right_index=True)
    
    count=0
    Only_ratio_n_mean_df=rmse_n_mean_n_std_n_ratio_df[["Observed Mean","Ratio_of_RMSE_and_Stdev"]]
    #Only_ratio_n_mean_df.plot.scatter("Observed Mean","Ratio_of_RMSE_and_Stdev")
    Only_ratio_n_mean_df.plot.bar()
    #for e in list(rmse_n_mean_n_std_n_ratio_df.index.values):
    #    mp.scatter(rmse_n_mean_n_std_n_ratio_df.loc(e,"Observed Mean"), rmse_n_mean_n_std_n_ratio_df.loc(e,"Ratio_of_RMSE_and_Stdev"), label=e, marker=count)
    #    count=count+1
    #    if(count==11):
    #        count=0
    rmse_n_mean_n_std_n_ratio_df.to_csv("/data/"+prefix+'_'+model_type+'_rmse_mean_std_and_ratio.csv')
    rmse_n_mean_n_std_n_ratio_n_r2_score_df.to_csv("/data/"+prefix+'_'+model_type+'_rmse_mean_std_and_ratio_and_r2_score.csv')
    #if len(list(rmse_n_mean_n_std_n_ratio_df.index.values)) <= 40:
    #    mp.legend(loc=(1.04,0))
    #else:
    #    mp.legend(bbox_to_anchor=(1.04, 1.04, 2.04, 2.04), loc='upper left', ncol=2, mode="expand")
    mp.xticks(rotation=90)
    mp.ylabel('Ratio_of_RMSE_and_Stdev')
    mp.xlabel(label_type+" Features")
    mp.title("Observed Mean and Ratio_of_RMSE_and_Stdev",size=12)
    mp.legend(loc=(1.04,0.5))
    #mp.savefig('test_result_plots_high_rmse.png',orientation='landscape',dpi=200,bbox_inches='tight')
    mp.savefig("/data/"+prefix+'_'+model_type+'_rmse_mean_std_and_ratio.png',bbox_inches='tight')
    mp.clf()

    Only_ratio_n_r2_score_df=rmse_n_mean_n_std_n_ratio_n_r2_score_df[["Ratio_of_RMSE_and_Stdev","r2_score"]]
    Only_ratio_n_r2_score_df.plot.bar()
    mp.xticks(rotation=90)
    mp.ylabel('r2_score')
    mp.xlabel(label_type+" Features")
    mp.title("Ratio_of_RMSE_and_Stdev and r2_score",size=12)
    mp.legend(loc=(1.04,0.5))
    #mp.savefig('test_result_plots_high_rmse.png',orientation='landscape',dpi=200,bbox_inches='tight')
    mp.savefig("/data/"+prefix+'_'+model_type+'_rmse_std_and_ratio_and_r2_score.png',bbox_inches='tight')
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
    #mp.legend(loc=(1.04,0))
    #mp.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
    #       ncol=2, mode="expand")
    if len(label_header_low_ratio) <= 40:
        mp.legend(loc=(1.04,0))
    else:
        mp.legend(bbox_to_anchor=(1.04, 1.04, 2.04, 2.04), loc='upper left', ncol=2, mode="expand")
    #mp.savefig('test_result_plots_low_rmse.png',orientation='landscape',dpi=200,bbox_inches='tight')
    mp.savefig("/data/"+prefix+'_'+model_type+'_test_result_plots_low_ratio.png',bbox_inches='tight')
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
    #mp.savefig('test_result_plots_high_rmse.png',orientation='landscape',dpi=200,bbox_inches='tight')
    mp.savefig("/data/"+prefix+'_'+model_type+'_test_result_plots_high_ratio.png',bbox_inches='tight')
    mp.clf()

    

    #Merging mean_df and rmse_dict
    if(len(ratio_low_dict)>1):
        ratio_low_df=pd.DataFrame.from_dict(ratio_low_dict,orient='index',columns=['RMSE/Stdev'])
        #mean_for_rmse_low_df=pd.DataFrame.from_dict(mean_for_rmse_low_dict,orient='index',columns=['Mean'])
        #mse_df=pd.DataFrame({'Label_features':label_list, 'MSE':mse_val})
        #print(mse_df)
        #mse_df.set_index('Label_features')
        ratio_low_df.to_csv("/data/"+prefix+"_"+model_type+"_Labels_with_Low_Ratio.csv")
        #mean_for_rmse_low_df.to_csv(prefix+"_"+model_type+"_Original_Mean_for_Labels_with_Low_RMSE.csv")

        #rmse_low_n_mean_df = rmse_low_df.merge(mean_for_rmse_low_df, how='outer', left_index=True, right_index=True)
        #mp.plot(rmse_low_df)
        ratio_low_df.plot(kind='bar')
        mp.ylabel('RMSE/Stdev')
        mp.xlabel('Label Features')
        mp.xticks(rotation=90)
        if len(label_header_low_ratio) <= 40:
            mp.xticks(size=5)
        else:
            mp.xticks(size=3)
        ##mp.title(i,size=1)
        mp.title("Low RMSE/Stdev for the label features",size=12)
        mp.savefig("/data/"+prefix+'_'+model_type+'_Low_Ratio_plot.png',orientation='landscape',dpi=100,bbox_inches='tight')
        mp.clf()

        #pd.plotting.scatter_matrix(rmse_low_n_mean_df)
        #mp.title("Comparing labels with Low RMSEs against their actual Mean values")
        #mp.savefig(prefix+'_'+model_type+'_Labels_w_Low_RMSE_and_actual_Means_plot.png',orientation='landscape',dpi=100,bbox_inches='tight')
        #mp.clf()

    else:
        print("Only 1 key:value pair in ratio_low_dict, so not proceeding with its plotting")
    if(len(ratio_high_dict)>1):
        ratio_high_df=pd.DataFrame.from_dict(ratio_high_dict,orient='index',columns=['RMSE/Stdev'])
        #mean_for_rmse_high_df=pd.DataFrame.from_dict(mean_for_rmse_high_dict,orient='index',columns=['RMSE'])
        #mse_df=pd.DataFrame({'Label_features':label_list, 'MSE':mse_val})
        #print(mse_df)
        #mse_df.set_index('Label_features')
        ratio_high_df.to_csv("/data/"+prefix+"_"+model_type+"_Labels_with_High_Ratio.csv")
        #mean_for_rmse_high_df.to_csv(prefix+"_"+model_type+"_Original_Mean_for_Labels_with_High_RMSE.csv")

        #rmse_high_n_mean_df = rmse_high_df.merge(mean_for_rmse_high_df, how='outer', left_index=True, right_index=True)

        mp.plot(ratio_high_df)
        mp.ylabel('RMSE/Stdev')
        mp.xlabel('Label Features')
        mp.xticks(rotation=90)
        mp.xticks(size=4)
        mp.title("High RMSE/Stdev for the label features",size=12)
        mp.savefig("/data/"+prefix+'_'+model_type+'_High_Ratio_plot.png',orientation='landscape',dpi=100,bbox_inches='tight')
        mp.clf()

        #pd.plotting.scatter_matrix(rmse_high_n_mean_df)
        #mp.title("Comparing labels with High RMSEs against their actual Mean values")
        #mp.savefig(prefix+'_'+model_type+'_Labels_w_High_RMSE_and_actual_Means_plot.png',orientation='landscape',dpi=100,bbox_inches='tight')
        #mp.clf()

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
    
    data_image1 = open("/data/"+prefix+'_'+model_type+'_test_result_plots_low_ratio.png', 'rb').read().encode('base64').replace('\n', '')
    img_tag1 = '<img src="data:image/png;base64,{0}">'.format(data_image1)
    outfileHTML.write(img_tag1+"\n")

    if(os.path.exists("/data/"+prefix+'_'+model_type+'_Low_Ratio_plot.png')):
        data_image2 = open("/data/"+prefix+'_'+model_type+'_Low_Ratio_plot.png', 'rb').read().encode('base64').replace('\n', '')
        img_tag2 = '<img src="data:image/png;base64,{0}">'.format(data_image2)
        outfileHTML.write(img_tag2+"\n")
    
    data_image3 = open("/data/"+prefix+'_'+model_type+'_rmse_mean_std_and_ratio.png', 'rb').read().encode('base64').replace('\n', '')
    img_tag3 = '<img src="data:image/png;base64,{0}">'.format(data_image3)
    outfileHTML.write(img_tag3+"\n")

    
    outfileHTML.write('Content-type: text/html\n\n'+""+"\n"+'<link href="default.css" rel="stylesheet" type="text/css" />')

    decision_thresholds=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

    ##Adding feature-importance calculations
    FI_handle=open("/data/"+prefix+"_"+"FeatureImportances.txt",'w+')
    X_features=list(test.columns)
    Y_features=list(Y_test.columns)
    if(model_type!="DecisionTree"):
        importance = model.coef_
        FI_handle.write("Y_Features, X_features_scored:"+str(X_features)+"\n")
        for u,k in enumerate(importance):
            FI_handle.write(str(Y_features[u])+","+str(k)+"\n")
    else:
    	importance = model.feature_importances_
        FI_handle.write("Features, Score"+"\n")
        for u,k in enumerate(importance):
            FI_handle.write(str(X_features[u])+","+str(k)+"\n")
    FI_handle.close()
    


    ##CALCULATE AUCs BASED on the devision thresholds list above.
    AUC_fH=open("/data/"+prefix+"_"+"AUC_values.txt",'w+')
    AUC_fH.write("label"+"\t"+"AUC_value"+"\t"+"decision_threshold"+"\n")
    label_resulted=[]
    for l in select_label_var_list:
        print l
        #fig, ax_=mp.subplots(1,1,figsize=(9,9))
        decision_threshold_list=[]
        AUC_value_list=[]
        #print Y_test_df[l]
        #print Y_test_df[[l]]
        for i in decision_thresholds:
            print i
            Y_test_df_n=binarize(pd.DataFrame.to_numpy(abs(Y_test_df[[l]])),threshold=i)
            Y_pred_df_n=binarize(pd.DataFrame.to_numpy(abs(Y_pred_df[[l]])),threshold=i)
            if(np.all((Y_test_df_n==0))):
                print("For decision threshold "+str(i)+":")
                print("Seems all values for label column "+l+" are zero. Hence not considering it for decision_threshold vs AUC plot")
            else:
                #fpr, tpr, thresholds = roc_curve(Y_test_df[[l]], Y_pred_df[[l]])
                fpr, tpr, thresholds = roc_curve(Y_test_df_n,Y_pred_df_n)
                #print thresholds
                #print fpr
                #print tpr
                #print auc(fpr,tpr)
                #print l; print i
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
                #if i==0.0:
                #    print "AUC for threshold "+str(i)+" and module "+l+"= "+str(AUC_value)
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
    #mp.xticks([0.0,0.2,0.4,0.6,0.8,1.0])
    mp.savefig("/data/"+prefix+'_'+model_type+'_AUC_for_decision_thresholds.png',orientation='landscape',dpi=100,bbox_inches='tight')
    mp.clf()
    data_image4 = open("/data/"+prefix+'_'+model_type+'_AUC_for_decision_thresholds.png', 'rb').read().encode('base64').replace('\n', '')
    img_tag4 = '<img src="data:image/png;base64,{0}">'.format(data_image4)
    outfileHTML.write(img_tag4+"\n")

    outfileHTML.write('Content-type: text/html\n\n'+""+"\n"+'<link href="default.css" rel="stylesheet" type="text/css" />')
    return label_resulted

## Model prediction function    
def predict(model , test):
    '''
    get model predictions
    '''
    Y_pred = model.predict(test)
    return Y_pred
##Defining the entire process starting from training through testing and validation per the modes of operation specified by user.       
def process(data_, label_, data_type, label_type, corr_method, corr_threshold, pVal_adjust_method, data_normalize_method, label_normalize_method, cv_par, scoring_par, mode, model_type, load_model, params, grid_search, param_grid, prediction_out, select_label_headers_for_predict, select_data_headers_for_predict):
    #outfileHTML=open(model_type+".output.html",'a')
    if os.path.isfile(data_):
        if os.path.getsize(data_)!=0:
            dataframe = read_dataset(data_)
            #print dataframe
        else:
            sys.exit("Size of data_ file"+data_+" is zero")
    else:
        sys.exit("data_ file:"+data_+" does not exists as a regular file")
    if mode != 'predict':
        if(os.path.isfile(label_)):
            if(os.path.getsize(label_)!=0):
                label = read_dataset(label_)
                #print label
                #label_header=list(label.keys())
            else:
                sys.exit("Size of label_ file:"+label_+"  is zero")
        else:
            sys.exit("data_ file:"+data_+" does not exists as a regular file")
    else:
        label='NA'
    #if(mode!='predict'):
    #    correlation(dataframe,label)
    dataframe , label, sampleIDs, label_header, dataframe_header =  preprocessing(dataframe , label, data_type, label_type, data_normalize_method, label_normalize_method, mode)
    #print "after dataframe"; print dataframe;
    #print "after label"; print label;
    
    if(mode!='predict' and mode!='validate'):
        select_data_var_list,select_label_var_list = correlation(dataframe,label,dataframe_header,label_header, model_type, corr_method, corr_threshold, pVal_adjust_method)
        print select_label_var_list
        #print label
        #print label[select_label_var_list]
        #print "filtered dataframe"; print filtered_dataframe
        #print(select_data_var_list)
        outfileHTML=open("/data/"+model_type+".output.html",'a')
        outfileHTML.write("<h1 style=text-align:center>"+"----------------------------------Features with highly significant correlations-------------------------------------"+"</h1>"+"\n")
        dataframe = dataframe[select_data_var_list]
        outfileHTML.write("<h2>"+"Below is the list of "+data_type+" features"+"</h2>"+"\n")
        outfileHTML.write("<h5>"+str(list(dataframe.keys()))+"</h5>"+"\n")
        #print("Below is the final list of data features used for training")
        #print(list(dataframe.keys()))
        label = label[select_label_var_list]
        print label

        outfileHTML.write("<h2>"+"Below is the list of "+label_type+" features"+"</h2>"+"\n")
        outfileHTML.write("<h5>"+str(list(label.keys()))+"</h5>"+"\n")
        #print("Below is the final list of label features used for training")
        #print(list(label.keys()))
        outfileHTML.write("<h2>"+"Number of "+data_type+" features"+"</h2>"+"\n")
        outfileHTML.write("<h3>"+str(len(select_data_var_list))+"</h3>"+"\n")
        #print("Below is the expected entries for data features")
        #print(len(select_label_var_list))
        outfileHTML.write("<h2>"+"Number of "+label_type+" features"+"</h2>"+"\n")
        outfileHTML.write("<h3>"+str(len(select_label_var_list))+"</h3>"+"\n")
        outfileHTML.close()
    elif(mode=='predict'):
        if(len(select_data_headers_for_predict)!=0):
            dataframe = dataframe[select_data_headers_for_predict]
        if(len(select_label_headers_for_predict)!=0):
            label_header_for_predict = select_label_headers_for_predict
    elif(mode=='validate'):
        if(len(select_data_headers_for_predict)!=0):
            dataframe = dataframe[select_data_headers_for_predict]
        if(len(select_label_headers_for_predict)!=0):
            select_label_var_list_for_validate=select_label_headers_for_predict
            label = label[select_label_var_list_for_validate]
        else:
            select_label_var_list_for_validate=label_header

    if mode == 'Train':
        train , Y_train ,test , Y_test = splitdata(dataframe , label, test_size, mode)
        print("Staring Training of :{}".format(model_type))
        model = BuildModel(train , Y_train , test , Y_test , model_type, params, cv_par, scoring_par, grid_search, param_grid, select_label_var_list, data_type, label_type, trainmodel = 'True')
        if save == 'True':
            try:
                joblib.dump(model, str(save_dir) + str(model_type)+ ".pkl" )
                print("Model saved at:{}".format(str(save_dir) + str(model_type)+ ".pkl"))
            except OSError:
                print("Saving model failed")
        else:
            print("Please provide save dir")
    
    elif mode == 'predict':
        print("Performing Prediction")
        try:
            model = joblib.load(str(load_model))
            #model = pickle.loads(load_model)
            print(model.get_params)

            outfileHTML=open("/data/"+model_type+".output.html",'a')
            store_params=model.get_params();
            outfileHTML.write("<h2 style=text-align:center;color:blue>"+"------------------------Model Summary-----------------------"+"</h2>")
            outfileHTML.write("<h3>"+"Model Parameters:"+"</h3>"+"\n")
            for i in store_params.keys():
                outfileHTML.write("<h4>"+str(i)+":"+str(store_params[i])+"</h4>")
            outfileHTML.write("<h2 style=text-align:center;color:green>"+"------------------------Samples for Prediction-----------------------"+"</h2>")
            outfileHTML.write("<h3>"+"No. of samples input for prediction: "+str(len(sampleIDs))+"</h3>")
            #outfileHTML.close()

            Y_pred = predict(model, dataframe)
            if prediction_out == "NULL":
               prediction_out = "prediction_out.txt"
            outfileH=open("/data/"+prediction_out,'w');

            outfileHTML.write("<h3>"+"All predicted values for labels available in the file: "+prediction_out+"</h3>")

            if(label_header_for_predict):
                outfileH.write("sampleID"+"\t"+"\t".join(label_header_for_predict)+"\n")
                outfileH.close()
                outfileH=open("/data/"+prediction_out,'a')
            #outfileHTML.write("<h4>"+"Checking if No. of samples provided as input = No. of samples predicted"+"</h4>")
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
    
    elif mode == 'validate':
        print("Performing validation")
        try:
            model = joblib.load(str(load_model))
            #model = pickle.loads(load_model)
            print(model.get_params)
            outfileHTML=open("/data/"+model_type+".output.html",'a')
            store_params=model.get_params();
            outfileHTML.write("<h2 style=text-align:center;color:blue>"+"------------------------Model Summary-----------------------"+"</h2>")
            outfileHTML.write("<h3>"+"Model Parameters:"+"</h3>"+"\n")
            for i in store_params.keys():
                outfileHTML.write("<h4>"+str(i)+":"+str(store_params[i])+"</h4>")
            outfileHTML.write("<h2 style=text-align:center;color:green>"+"------------------------Samples for Validation-----------------------"+"</h2>")
            outfileHTML.write("<h3>"+"No. of samples used for validation: "+str(len(sampleIDs))+"</h3>")
            outfileHTML.close()
            #for i in label_header:
            #dataframe = dataframe.loc[:, (dataframe<=0).any(axis=0)]
            evaluate(model, dataframe , label, select_label_var_list_for_validate, 'validation', model_type, data_type, label_type)
        except IOError:
            print("Probably saved model file is not provided/check path")
    else:
        sys.exit("Invalid mode. Should be either Train, predict or validate. Please check your config file.")

## MAIN function. The main starting point of the script. Scans parameters, defines variable and passes on to the process function called below.
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
    #classification = config['modes']['classification']
    ##getting select data and label features from respective parameters in 'modes' section
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
    #print select_data_headers_for_predict; print select_label_headers_for_predict
    
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
    outfileHTML.write("<h1 style=text-align:center;color:red;>"+"Radiogenomics Analysis Report"+"</h1>"+"\n")
    #write date and time of the report
    from datetime import datetime
    datetime_now = datetime.now()
    datetime_string = datetime_now.strftime("%d/%m/%Y %H:%M:%S")
    outfileHTML.write("<h5 style=text-align:center;>"+datetime_string+"</h5>"+"\n")
    outfileHTML.write("<h2 style=text-align:center;color:brown>"+"----------------------------Model inputs-------------------------------"+"</h2>")
    outfileHTML.write("<h3>"+"Mode:{}".format(mode)+"</h3>")
    outfileHTML.write("<h3>"+"Model:{}".format(model_type)+"</h3>")
    outfileHTML.write("<h3>"+"Params:{}".format(print_params)+"</h3>")
    if grid_search == "True":
        outfileHTML.write("<h3>"+"Grid_Params:{}".format(print_param_grid)+"</h3>")
    outfileHTML.close()
    
    ## CALLING THE PROCESS FUNCTION.


    process(args.data, args.label, data_type, label_type, corr_method, corr_threshold, pVal_adjust_method, data_normalize_method, label_normalize_method, cv_par, scoring_par, mode, model_type, load_model, params, grid_search, param_grid, args.prediction_out, select_label_headers_for_predict, select_data_headers_for_predict)
