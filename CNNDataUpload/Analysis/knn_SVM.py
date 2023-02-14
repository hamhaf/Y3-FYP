#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 00:47:28 2020

@author: shariba
"""

import pandas as pd
import numpy as np
import datatable as dt
import time
from sklearn.impute import SimpleImputer
from utils.metric_computation import test_predict
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
import os
import argparse

debug = 0

# plots the confusion matrices
def plots_prints_metrics(clf_svm, y_test1, predicted, resultFileName, nclass):
    if nclass==2:
        from utils.plot_PR_confusionMat import confusionPlot  
        cfMat = confusionPlot (y_test1, predicted, resultFileName, str(classA), str(classB))
    else:
        from utils.plot_PR_confusionMat_3class import confusionPlot 
        cfMat = confusionPlot (y_test1, predicted, resultFileName)
    print("classification result for %s:\n%s\n"
          % (clf_svm, metrics.classification_report(y_test1, predicted)))
    print("Confusion matrix:\n%s" % metrics.confusion_matrix(y_test1, predicted))
    return cfMat

# not entirely sure about this normalisation function
def normaliseLabels(ylabel, nlabels):

    if (nlabels == 3):
         ylabel = ylabel-1
    else:
         if (ylabel.max() == 2):
             ylabel = ylabel-1
         else:
             for idx in range(len(ylabel)):
                 if (ylabel[idx][0] == 3):
                     ylabel[idx] = ylabel[idx]-2 
                 elif (ylabel[idx][0] == 2):
                     ylabel[idx] = ylabel[idx]-2 
                 else:
                     ylabel[idx] = ylabel[idx]-1
    return ylabel
             
             
# gets the arguments from cmd line
def get_argparser():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_data", type=str, default='./data/noExclusion_train_data.csv',
                        help="path to train")
    parser.add_argument("--train_label", type=str, default='./data/noExclusion_train_label.csv',
                        help="path to train")
    
    parser.add_argument("--test_data", type=str, default='./data/noExclusion_test_data.csv',
                        help="path to train")
    parser.add_argument("--test_label", type=str, default='./data/noExclusion_test_label.csv',
                        help="path to train")
    
    parser.add_argument("--num_labels", type=int, default=3,
                        help="num classes (default: 2/3 class)")
    
    parser.add_argument("--outputFolder_ckpt", type=str, default='./results_Classical_ML/',
                        help="output folder")
    
    parser.add_argument("--gpu_id", type=str, default='0',
                        help="GPU ID")
    
    parser.add_argument("--classSplit", type=str, default='2,3',
                        help="1,2")
    
    return parser



if __name__ == '__main__':
    
    opts = get_argparser().parse_args()
    # gets filename from path by taking last element of 1st split, then removes file extension
    # from filename with second split
    fileName = ((opts.train_data).split('/')[-1]).split('.')[0]
    
    classA = int(opts.classSplit.split(',')[0])
    classB = int(opts.classSplit.split(',')[1])
    
    val = 0
    nlabel=opts.num_labels

    """prepare train and test data"""
    csv_test_data =  opts.test_data
    csv_test_labels  =  opts.test_label
        
    csv_train_data =  opts.train_data
    csv_train_labels  =  opts.train_label
    if nlabel==2:
        from utils.plot_PR_confusionMat import confusionPlot  
    else:
        from utils.plot_PR_confusionMat_3class import confusionPlot 

        
    # make result folder
    directoryName = opts.outputFolder_ckpt
    os.makedirs(directoryName, exist_ok = True) 
    
    # data loading 
    data = pd.read_csv(csv_train_data, header = None)
    df = pd.DataFrame(data.values)
    # convert df to np array - X_train1 will be an array of arrays (2D array)
    # where each element represents a row from the dataset
    X_train1 = np.array(df, dtype=float)
    data = pd.read_csv(csv_test_data, header = None)
    df = pd.DataFrame(data.values)
    X_test1 = np.array(df, dtype=float)
    # read in the training labels
    ytrain = pd.read_csv(csv_train_labels, header=None)
    # convert to 1D np array
    y_en= np.array(ytrain[0])
    # reshape into 1D normal array
    y_train1 = y_en.reshape(y_en.shape[0],1) 
    # normalise labels
    y_train1 = normaliseLabels(y_train1, nlabel)
    # for two class we need to adject no.

    y_test = pd.read_csv(csv_test_labels, header=None)
    y_en= np.array(y_test[0])
    y_test1 = y_en.reshape(y_en.shape[0],1)
    
    y_test1 = normaliseLabels(y_test1, nlabel)

    # perform PCA
    # what is PCA: dimensionality-reduction method that is often used to 
    # reduce the dimensionality of large data sets, by transforming a large set 
    # of variables into a smaller one that still contains most of the information in
    # the large set.
    from sklearn.decomposition import PCA
    time_start = time.time()
    # Ncomponents=300 # number of coordinates
    #your code here (store the components in X_reduced)
    # pca = PCA(n_components=Ncomponents)
    pca = PCA()
    X_reduced=pca.fit_transform(X_train1)
    print(X_reduced.shape)
    print('PCA time: {} seconds'.format(time.time()-time_start))
    
    ## LDA/QDA
    from sklearn.model_selection import train_test_split
    # X_train1, X_test1, y_train1, y_test1 = train_test_split(X_norm, y, test_size = 0.20, random_state =0)
    print("Training set with feature selection has {} samples.".format(X_train1.shape[0]))
    print("Testing set with feature selection has {} samples.".format(X_test1.shape[0]))
    
    
    if nlabel==3:
        fileName = fileName.split('_')[0]+'_3way'
        jsonFileName = fileName+ '_3way.json'

    else:
        fileName = fileName.split('_')[0] +'_classA_'+str(classA)+'_classB_'+str(classB)+'_2way'
        jsonFileName = fileName+'_2way.json'

    
    #========================= KNN ============================
    trainScores={}
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.model_selection import cross_val_score
    bestScore=0.0
    accList=[]
    for k in range(1,10):
        clf_knn = KNeighborsClassifier(n_neighbors=k,algorithm='auto')
        # using 10 fold cross validation for scoring the classifier's accuracy
        scores = cross_val_score(clf_knn, X_train1, y_train1, cv=5)
        score=scores.mean()
        print(score)
        accList.append(score)
        if score > bestScore:
            bestScore=score
            best_clf=clf_knn
            bestK=k
    print("Best K is :",bestK,"| Cross validation Accuracy :",bestScore)
    clf_knn=best_clf
    ## validate with your test samples
    if debug:
        plt.plot(range(1,10),accList)
        plt.xlabel('K')
        plt.ylabel('CV Accuracy')
        plt.show()
    # your code here
    clf_knn.fit(X_train1,y_train1)
    y_pred=best_clf.predict(X_train1)
    # find train accuracy
    print('Accuracy of KNeighborsClassifier on training set: {:.2f}'.format(best_clf.score(X_train1, y_train1)))
    val_KNN=test_predict(best_clf, X_test1, y_test1)
    # Confusion matrixs
    predicted = best_clf.predict(X_test1)
    # print("classification result for %s:\n%s\n"
    #       % (best_clf, metrics.classification_report(y_test1, predicted)))
    # print("Confusion matrix:\n%s" % metrics.confusion_matrix(y_test1, predicted))
    
    cfMat_KNN = plots_prints_metrics(best_clf, y_test1, predicted, directoryName+'/'+fileName + '_KNN_best_k'+ str(bestK) + '.svg', nlabel)

    #========================= SVM ==================================
    from sklearn import svm
    # ply works better than rbf and linear
    clf_svm = svm.SVC(kernel='rbf', gamma='scale', C=1000)
    clf_svm.fit(X_train1, y_train1.ravel())
    y_pred = clf_svm.predict(X_test1)
    print('Accuracy of SVM on training set: {:.2f}'.format(clf_svm.score(X_train1, y_train1)))
    val_SVM=test_predict(clf_svm, X_test1, y_test1)
    predicted = clf_svm.predict(X_test1)
    #==================================================================================
    #==================================================================================
    y_score = clf_svm.decision_function(X_test1)
    # plot_PR_curve(y_test1, y_score, predicted, directoryName+'/'+fileName + '_SVM_PR.svg')
    
    cfMat_svm = plots_prints_metrics(clf_svm, y_test1, predicted, directoryName+'/'+fileName + '_SVM.svg', nlabel)
    # svm classwise
    from utils.multi_classwise_metric import compute_cfMatrix, compute_metricswithCFMatrix, saveMetrics2Json
    
    sensitivity, specificity, PPV, NPV, Accuracy, F1 = compute_cfMatrix(predicted, y_test1, labels=[1, 2, 3])
    # labels=[1, 2, 3, 4]
    val_SVM=test_predict(clf_svm, X_test1, y_test1)
    print ('Accuracy per class average:', np.mean(Accuracy))
    print ('F1 per class average:', np.mean(F1))

    """
    Evaluate metrics and write it on json file
    """
    import os
    import json
    Result_dir=directoryName
    os.makedirs(Result_dir, exist_ok=True)
    jsonFileName=os.path.join(Result_dir, jsonFileName)

    
    sensitivity, specificity, PPV, NPV, Accuracy,FPR, FNR, FDR, F1 = compute_metricswithCFMatrix(cfMat_KNN, nlabel)
    my_dictionary = saveMetrics2Json('KNN',sensitivity, specificity, PPV, NPV, Accuracy,FPR, FNR, FDR, F1)
    fileObj= open(jsonFileName, "a")
    fileObj.write("\n")
    json.dump(my_dictionary, fileObj)
    ############################
    fileObj.write("\n")
    sensitivity, specificity, PPV, NPV, Accuracy,FPR, FNR, FDR, F1 = compute_metricswithCFMatrix(cfMat_svm, nlabel)
    my_dictionary = saveMetrics2Json('SVM',sensitivity, specificity, PPV, NPV, Accuracy,FPR, FNR, FDR, F1)
    json.dump(my_dictionary, fileObj)
    
    fileObj.close()
    print('Evaluation metric values written to---->', (jsonFileName))
    

    # data statistics
    
