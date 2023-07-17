#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 12:07:53 2020
@author: sharib
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SubsetRandomSampler
import argparse

from model.model import model

        
def get_argparser():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_data", type=str, default='./data/noExclusion_train_data.csv',
                        help="path to train")
    parser.add_argument("--train_label", type=str, default='./data/noExclusion_train_label.csv',
                        help="path to train")
    parser.add_argument("--num_labels", type=int, default=3,
                        help="num classes (default: 2/3 class)")
    
    parser.add_argument("--batch_size", type=int, default=16,
                        help='batch size (default: 16)')
    parser.add_argument("--val_batch_size", type=int, default=8,
                        help='batch size for validation (default: 4)')
    parser.add_argument("--trainsplit", type=float, default=0.20,
                        help="train-val split")
    
    parser.add_argument("--outputFolder_ckpt", type=str, default='./ckpt/',
                        help="output folder")
    
    parser.add_argument("--gpu_id", type=str, default='0',
                        help="GPU ID")
                        
    parser.add_argument("--classSplit", type=str, default='2,3',
                                            help="1,2")
    
    return parser



if __name__ == '__main__':
    import numpy as np
    import torch.optim as optim
    import os
    import time
    from sklearn.model_selection import KFold, StratifiedKFold
    
    
    opts = get_argparser().parse_args()
    
    fileName = ((opts.train_data).split('/')[-1]).split('.')[0]
    
    os.makedirs(opts.outputFolder_ckpt, exist_ok= True)
    
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    # print(device)
    # if device !='cpu':
        # print('cuda detected -------->')
        # os.environ["CUDA_VISIBLE_DEVICES"] =opts.gpu_id
        # device = 'cuda'
        # torch.backends.cudnn.deterministic = True
        # torch.backends.cudnn.benchmark = False
    print(device)
    nlabel = opts.num_labels

    classSplit = opts.classSplit
    classA = int(classSplit.split(',')[0])
    classB = int(classSplit.split(',')[1])
    
    # seeds placed for reproducibility
    seed = 0
    np.random.seed(seed)


    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(0)
    
    m = model(nlabel).to(device)

    # best parameter 
    lr = 0.001
    decay_rate = 0.0001
    beta_1 = 0.92
    beta_2 = 0.999
    
    if fileName.split('_')[0] == 'balanced':
        print('balanced data detected...using tuned hyperparameters...')
        lr = 0.0008
        if nlabel == 3:
            decay_rate = 0.00001


    optimizer = optim.Adam(m.parameters(), betas=(beta_1, beta_2), lr=lr, weight_decay=decay_rate*lr)
    criterion = nn.NLLLoss()
    epochs = 500
    if nlabel==3:
        from utils.dataLoader_HSI_sept_3class import MUSE_HSI_dataLoader
        hsi_dataset = MUSE_HSI_dataLoader(csv_file_data=opts.train_data, csv_file_labels=opts.train_label)
    else:
        from utils.dataLoader_HSI_2class import MUSE_HSI_dataLoader
        hsi_dataset = MUSE_HSI_dataLoader(csv_file_data=opts.train_data, csv_file_labels=opts.train_label, class_categories=classSplit)
   
    useValLoss = 1

    # scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'max')
    
    if useValLoss == 1:
        validation_split = opts.trainsplit
        percentageSplit = validation_split * 100
        shuffle_dataset = False
        dataset_size = len(hsi_dataset)
        indices = list(range(dataset_size))
        split = int(np.floor(validation_split * dataset_size))
        
        if shuffle_dataset :
            np.random.seed(seed)
            # np.random.shuffle(indices)
            
        train_indices, val_indices = indices[split:], indices[:split]
            
        # # Creating PT data samplers and loaders:
        train_sampler = SubsetRandomSampler(train_indices)
        valid_sampler = SubsetRandomSampler(val_indices)
    
        print('length of train data:', len(train_sampler))
        print('length of val data:', len(valid_sampler))
        # here batch 16 is the right choice
        dataloader_train = DataLoader(hsi_dataset, batch_size=opts.batch_size, sampler=train_sampler,  num_workers=0) 
        dataloader_val = DataLoader(hsi_dataset, batch_size=opts.val_batch_size, sampler=valid_sampler, num_workers=0) 
        
    else:
        kf = KFold(n_splits = 5)
        skf = StratifiedKFold(n_split = 5, random_state = 7, shuffle = True) 
        hsi_dataset = MUSE_HSI_dataLoader(csv_file_data=opts.train_data, csv_file_labels=opts.train_label)
        dataloader_train = DataLoader(hsi_dataset, batch_size=opts.batch_size,   shuffle=True, num_workers=1) 

    print_every = 10
    test_losses = []
    train_losses = []  
    running_loss = 0

    best_acc = 0
    best_loss = 100.0
    current_loss = 0.0
    st = time.time()
    for epoch in range(epochs):
        losses = []
        test_loss = 0
        accuracy_train = 0
        for i, sample in enumerate(dataloader_train):
            if device=='cuda':
                inputv = torch.cuda.FloatTensor(sample['feature'].to(device))
                labelsv = torch.cuda.FloatTensor(sample['label'].to(device).type(torch.float).unsqueeze(1))
            else:
                inputv = torch.FloatTensor(sample['feature'].to(device))
                labelsv = torch.FloatTensor(sample['label'].to(device).type(torch.float).unsqueeze(1))
            
            output = m(inputv)
            loss = criterion(output, labelsv.type(torch.long).flatten())
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            running_loss += loss.data.item()
            losses.append(loss.data.item())
          
        if useValLoss == 1:
            
            if epoch%print_every == 0:
                accuracy = 0
                m.eval()
                with torch.no_grad():
                    for k, sample_val in enumerate(dataloader_val):
                        if device=='cuda':
                            inputs = torch.cuda.FloatTensor(sample['feature'].to(device))
                            labels = torch.cuda.FloatTensor(sample['label'].to(device).type(torch.float).unsqueeze(1))
                        else:
                            inputs = torch.FloatTensor(sample['feature'].to(device))
                            labels = torch.FloatTensor(sample['label'].to(device).type(torch.float).unsqueeze(1))
                        logps = m.predict(inputs)
                        
                        batch_loss = criterion(logps, labels.type(torch.long).flatten())
                        test_loss += batch_loss.item()
                        ps = torch.exp(logps)
                        top_p, top_class = ps.topk(1, dim=1)
    
                        equals = top_class == labels.view(*top_class.shape)
                        accuracy +=torch.mean(equals.type(torch.FloatTensor)).item()
                        
                        # scheduler.step(batch_loss.item())
                       
                train_losses.append(running_loss/len(dataloader_train))
                test_losses.append(test_loss/len(dataloader_val))
                
                print(f"Epoch {epoch+1}/{epochs}.. "
                      f"Train loss: {running_loss/print_every:.6f}.. "
                      # f"Train accuracy: {accuracy_train/print_every:.6f}.."
                      f"Validation loss: {test_loss/len(dataloader_val):.6f}.. "
                      f"Validation accuracy: {accuracy/len(dataloader_val):.6f}")
               
                acc = accuracy/len(dataloader_val)
                current_loss = float(running_loss/print_every)
                if (acc >= best_acc) and (epoch>200) and (current_loss < best_loss): #and acc  <1.0
                    best_acc = acc
                    print('saving ckpt...')
                    # 2-way classification stores with classs categories unlike 3-way
                    if nlabel==2:
                        torch.save(m.state_dict(), os.path.join(opts.outputFolder_ckpt, fileName.split('_')[0] + '_sept_nclasses_%d_v%d_%d_%d.pth'%( (nlabel),  (percentageSplit), classA, classB )))
                    else:
                        torch.save(m.state_dict(), os.path.join(opts.outputFolder_ckpt, fileName + '_sept_nclasses_%d_v%d.pth'%( (nlabel),  (percentageSplit) )))
                    
                if  (current_loss < best_loss):
                    best_loss = current_loss
                running_loss = 0

    # save the final checkpoint
    if nlabel==3:
        torch.save(m.state_dict(), os.path.join(opts.outputFolder_ckpt, fileName + 'nclasses_final_%d_v%d.pth'%( (nlabel), (percentageSplit) )))
    else:
        torch.save(m.state_dict(), os.path.join(opts.outputFolder_ckpt, fileName.split('_')[0] + 'nclasses_final_%d_v%d_%d_%d.pth'%( (nlabel), (percentageSplit), classA, classB )))
        
    print('Elapsed training time: {}'.format(time.time()-st))
    
    debug = 0
#    If debug then plot train-val loss
    if debug:
        import matplotlib.pyplot as plt
        plt.plot(train_losses, label='Training loss')
        plt.plot(test_losses, label='Validation loss')
        plt.legend(frameon=False)
        plt.show()
