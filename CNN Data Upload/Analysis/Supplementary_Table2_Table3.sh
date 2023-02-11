#!/bin/bash
# Classical Machine Learning Methods: knn and svm methods (Table 2 and Table 3)
# For both 3 class classification and 2 class classification
# Results will be written into a new folder Results_
# input: training data and labels
# output: confusion matrix written in folder

BASE_FOLDER=`pwd`
echo $BASE_FOLDER
ROOT_FOLDER=$BASE_FOLDER

# Step1: conda environment (set your environment)
    #  check how to install your virtual environment here: https://github.com/sharibox/tutorial/blob/master/conda_essentials.md
    # Install dependencies: scikit-learn, matplotlib

# Step 2: Activate your environment
source activate myenv

# Step3: Set class and your dataset
nClasses=3
dataType='noExclusion' 

if (( $nClasses == '3' ))
then 
##========================= Classical ML for 3 class category training ===============================================##
    echo '---->' $dataType 'data for training... for ' $nClasses'-way classification'
    python $BASE_FOLDER/knn_SVM.py --train_data $BASE_FOLDER/data/$dataType'_train_data.csv' --train_label $BASE_FOLDER/data/$dataType'_train_label.csv' \
    --test_data $BASE_FOLDER/data/$dataType'_test_data.csv' --test_label $BASE_FOLDER/data/$dataType'_test_label.csv' --num_labels 3  --classSplit '1,2,3'

fi 

nClasses=2
##========================= Classical ML 2 class category training ===============================================##
if (( $nClasses == 2 ))
then 

    echo '---->' $dataType 'data for training... for ' $nClasses'-way classification'
     
    classPairSrting='1,2'
    classA='1'
    classB='2'
    python $BASE_FOLDER/knn_SVM.py --train_data $BASE_FOLDER/data/$dataType'_train_data_'$classA'_'$classB'.csv' --train_label $BASE_FOLDER/data/$dataType'_train_label_'$classA'_'$classB'.csv' \
    --test_data $BASE_FOLDER/data/$dataType'_test_data_'$classA'_'$classB'.csv' --test_label $BASE_FOLDER/data/$dataType'_test_label_'$classA'_'$classB'.csv' \
     --num_labels 2  --classSplit $classPairSrting

    classPairSrting='1,3'
    classA='1'
    classB='3'
    python $BASE_FOLDER/knn_SVM.py --train_data $BASE_FOLDER/data/$dataType'_train_data_'$classA'_'$classB'.csv' --train_label $BASE_FOLDER/data/$dataType'_train_label_'$classA'_'$classB'.csv' \
    --test_data $BASE_FOLDER/data/$dataType'_test_data_'$classA'_'$classB'.csv' --test_label $BASE_FOLDER/data/$dataType'_test_label_'$classA'_'$classB'.csv' \
     --num_labels 2  --classSplit $classPairSrting


    classPairSrting='2,3'
    classA='2'
    classB='3'
    python $BASE_FOLDER/knn_SVM.py --train_data $BASE_FOLDER/data/$dataType'_train_data_'$classA'_'$classB'.csv' --train_label $BASE_FOLDER/data/$dataType'_train_label_'$classA'_'$classB'.csv' \
    --test_data $BASE_FOLDER/data/$dataType'_test_data_'$classA'_'$classB'.csv' --test_label $BASE_FOLDER/data/$dataType'_test_label_'$classA'_'$classB'.csv' \
     --num_labels 2  --classSplit $classPairSrting

fi
