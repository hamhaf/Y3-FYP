#!/bin/bash

BASE_FOLDER=`pwd`
echo $BASE_FOLDER
ROOT_FOLDER=$BASE_FOLDER
source activate myenv

# for 3 class training use 3 else use 2
nClasses=3
dataType=0 # selection 0 , 1

#echo '================> Training begins for 3-way classification <================================'
if (( $nClasses == '3' ))
then 
##========================= 3 class category training ===============================================##
# TO run 3 class classification using CNN on "No exclusion and balanced" dataset

   echo '----> noExclusion data for training...'
   python $BASE_FOLDER/train.py --train_data $BASE_FOLDER/data/noExclusion_train_data.csv --train_label $BASE_FOLDER/data/noExclusion_train_label.csv \
   --num_labels $nClasses --batch_size 16 --val_batch_size 8 --trainsplit 0.20 --outputFolder_ckpt $BASE_FOLDER/ckptDir --gpu_id 0
   
    echo '----> balanced data for training...'
   # TO run 3 class classification using CNN on "balanced" dataset
   python $BASE_FOLDER/train.py --train_data $BASE_FOLDER/data/balanced_train_data.csv --train_label $BASE_FOLDER/data/balanced_train_label.csv \
   --num_labels $nClasses --batch_size 16 --val_batch_size 8 --trainsplit 0.20 --outputFolder_ckpt $BASE_FOLDER/ckptDir --gpu_id 0

fi 

nClasses=2
#echo '================> Training begins for 2-way classification <================================'
###========================= 2 class category training ===============================================##
if (( $nClasses == 2 ))
then

       echo '----> noExclusion data for training for class 1 and 2...'
       python $BASE_FOLDER/train.py --train_data $BASE_FOLDER/data/noExclusion_train_data_1_2.csv --train_label $BASE_FOLDER/data/noExclusion_train_label_1_2.csv \
       --num_labels $nClasses --batch_size 16 --val_batch_size 4 --trainsplit 0.20 --outputFolder_ckpt $BASE_FOLDER/ckptDir --gpu_id 0 --classSplit '1,2'


       echo '----> noExclusion data for training for class 1 and 3...'
       python $BASE_FOLDER/train.py --train_data $BASE_FOLDER/data/noExclusion_train_data_1_3.csv --train_label $BASE_FOLDER/data/noExclusion_train_label_1_3.csv \
       --num_labels $nClasses --batch_size 16 --val_batch_size 4 --trainsplit 0.20 --outputFolder_ckpt $BASE_FOLDER/ckptDir --gpu_id 0 --classSplit '1,3'


       echo '----> noExclusion data for training for class 2 and 3...'
       python $BASE_FOLDER/train.py --train_data $BASE_FOLDER/data/noExclusion_train_data_2_3.csv --train_label $BASE_FOLDER/data/noExclusion_train_label_2_3.csv \
       --num_labels $nClasses --batch_size 16 --val_batch_size 4 --trainsplit 0.20 --outputFolder_ckpt $BASE_FOLDER/ckptDir --gpu_id 0 --classSplit '2,3'


    #    else for balanced!
         echo '----> balanced data for training for class 1 and 2...'
        python $BASE_FOLDER/train.py --train_data $BASE_FOLDER/data/balanced_train_data_1_2.csv --train_label $BASE_FOLDER/data/balanced_train_label_1_2.csv \
        --num_labels $nClasses --batch_size 16 --val_batch_size 4 --trainsplit 0.20 --outputFolder_ckpt $BASE_FOLDER/ckptDir --gpu_id 0 --classSplit '1,2'

        echo '----> balanced data for training for class 1 and 3...'

        python $BASE_FOLDER/train.py --train_data $BASE_FOLDER/data/balanced_train_data_1_3.csv --train_label $BASE_FOLDER/data/balanced_train_label_1_3.csv \
        --num_labels $nClasses --batch_size 16 --val_batch_size 4 --trainsplit 0.20 --outputFolder_ckpt $BASE_FOLDER/ckptDir --gpu_id 0 --classSplit '1,3'

        echo '----> balanced data for training for class 2 and 3...'
        python $BASE_FOLDER/train.py --train_data $BASE_FOLDER/data/balanced_train_data_2_3.csv --train_label $BASE_FOLDER/data/balanced_train_label_2_3.csv \
        --num_labels $nClasses --batch_size 16 --val_batch_size 4 --trainsplit 0.20 --outputFolder_ckpt $BASE_FOLDER/ckptDir --gpu_id 0 --classSplit  '2,3'

fi
echo '================> End of all training, you can now run test script to reproduce results<======================='
