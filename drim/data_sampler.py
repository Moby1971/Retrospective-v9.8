__author__ = 'BMEP'

import os
import numpy as np
from torch.utils.data import Dataset
import logging
from scipy.io import loadmat
import h5py
import matplotlib.pyplot as plt
import numpy.fft as fft
import imageio
import cv2
import math


logger = logging.getLogger(__name__)

class MRData(Dataset):

    def __init__(self, data_path):

        self.data_path = data_path
        self.subjects = ['retroAItemp.mat']    
        #self.subjects = ['retroData_1_12_24043_000_0_f54.mat']
        logger.info('Subject included:')

        self.crop = 100
        
        lengths = []
        self.data = dict()
        for subject in self.subjects:
            print(subject)
    
            logger.info(subject)
            self.data[subject] = dict()
            
            mat_contents = loadmat(self.data_path + '/' + subject)
            kspace = np.array(mat_contents['kData'])
            undersampled_kspace_mask = np.array(mat_contents['fData'])
        
            if len(kspace.shape) == 4:
                kspace = kspace[0, :, :, :]
                undersampled_kspace_mask = np.expand_dims(undersampled_kspace_mask, axis=0)
                kspace = np.expand_dims(kspace, axis=0)

            if len(kspace.shape) == 5:
                kspace = kspace[0, :, :, :, :]
                undersampled_kspace_mask = undersampled_kspace_mask.transpose(3, 0, 2, 1)
                kspace = kspace.transpose(3, 0, 1, 2)
      
            print(f'Undersampled kspace mask shape: {undersampled_kspace_mask.shape}')
            print(f'Kspace shape: {kspace.shape}')
         
            self.data[subject]['kspace'] = kspace
            self.data[subject]['mask'] = undersampled_kspace_mask
            self.data[subject]['sense'] = np.ones_like(kspace)
            
            lengths.append(kspace.shape[2])

        logger.info('Finished pre-processing subject ...')

        self.lengths = np.cumsum(lengths)



    def __len__(self):
        return self.lengths[-1]



    def __getitem__(self, index):

        idx = np.digitize(index, self.lengths)
        if idx > 0:
            index -= self.lengths[idx - 1]
        subject = self.subjects[idx]

        imspace = fft.ifftshift(fft.ifft2(self.data[subject]['kspace'],), axes = (2,3))
        imspace = cv2.normalize(np.abs(imspace), None, 0, 1, cv2.NORM_MINMAX)
        kspace_undersampled_mask = self.data[subject]['mask']
        sense = self.data[subject]['sense'] # [C, D, H, W]

        imspace = imspace[:, :, index] # [D, dyn, W]
        kspace_undersampled_mask = kspace_undersampled_mask[:, :, index]    
        sense = sense[:, :, index] # [C, D, W]
            
        imspace = imspace.transpose(1, 0, 2) # [T, D, W]
        kspace_undersampled_mask = kspace_undersampled_mask.transpose(1, 0, 2)
        sense = sense.transpose(1, 0, 2)
       
        # DATA 
        sense = sense[np.newaxis] # [C, T, D, W]
        imspace = imspace[np.newaxis]
        kspace_undersampled_mask = kspace_undersampled_mask.astype(np.int8)
        kspace = np.fft.fft(sense * imspace, axis=-1)
        estimate = np.sum(np.fft.ifft(kspace, axis=-1) * sense.conj(), 0)

        # print(estimate[0,0,:])

        #print(f'sense shape: {sense.shape}')
        #print(f'mask shape: {kspace_undersampled_mask[np.newaxis, ..., np.newaxis].shape}')
        #print(f'measurement shape: {kspace.shape}')
        #print(f'estimate shape: {kspace.shape}')

        data = {
            'sense': np.ascontiguousarray(sense), # [C, T, D, W]
            'mask': kspace_undersampled_mask[np.newaxis, ..., np.newaxis], # [1, T, D, W, 1]
            'measurements': kspace, # [C, T, D, W]
            'estimate': estimate # [T, D, W]
        }

        return data
