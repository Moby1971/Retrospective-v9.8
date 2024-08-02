__author__ = 'Kai'

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


logger = logging.getLogger(__name__)

class MRData(Dataset):
    def __init__(self, data_path, crop, train):

        self.train = train
        self.data_path = data_path
        self.subjects = ['retroAItemp.mat']    # os.listdir(self.data_path)
        logger.info('Subject included in training:')

        self.crop = crop
        
        lengths = []
        self.data = dict()
        for subject in self.subjects:
            print(subject)
    
            logger.info(subject)
            self.data[subject] = dict()
            
            try:
                mat_contents = loadmat(self.data_path + '/' + subject)
                kspace = np.array(mat_contents['kData'])
                undersampled_kspace_mask = np.array(mat_contents['fData'])
                    
                if len(kspace.shape) == 4:
                    kspace = kspace[0, :, :, :]
                    undersampled_kspace_mask = np.expand_dims(undersampled_kspace_mask, axis=0)
                    kspace = np.expand_dims(kspace, axis=0)

                if len(kspace.shape) == 5:
                    kspace = kspace[0, :, :, :, :]
                    undersampled_kspace_mask = undersampled_kspace_mask.transpose(3, 0, 1, 2)
                    kspace = kspace.transpose(3, 0, 1, 2)

                print(f'Undersampled kspace mask shape: {undersampled_kspace_mask.shape}')
                print(f'Kspace shape: {kspace.shape}')
                print('Matlab')
        
            except Exception as e:
                with h5py.File(self.data_path + '/' + subject, 'r') as file:
                    kspace = np.array(file['kData'])
                    undersampled_kspace_mask = np.array(file['fData'])
                kspace = kspace.transpose(3, 2, 1, 0)
                real = kspace['real']
                imag = kspace['imag']
                kspace = real + 1j * imag
                undersampled_kspace_mask = undersampled_kspace_mask[0, :, :, :, :]
                undersampled_kspace_mask = undersampled_kspace_mask.transpose(0, 3, 2, 1)
                print(undersampled_kspace_mask.shape)
                print(f'Kspace shape: {kspace.shape}')

            self.data[subject]['kspace'] = kspace
            self.data[subject]['mask'] = undersampled_kspace_mask
            self.data[subject]['sense'] = np.ones_like(kspace)

            
            lengths.append(kspace.shape[2])

        logger.info(
            f'Finished pre-processing subject{"s" if self.train else ""}.')

        self.lengths = np.cumsum(lengths)


    def __len__(self):
        return self.lengths[-1]

    def __getitem__(self, index):
        idx = np.digitize(index, self.lengths)
        if idx > 0:
            index -= self.lengths[idx - 1]
        subject = self.subjects[idx]
        kspace = self.data[subject]['kspace'] # [D, dyn, H, W]
        imspace = fft.ifftshift(fft.ifft2(self.data[subject]['kspace'],), axes = (2,3)) # axes 2,3
        
        imspace_undersampled_mask = self.data[subject]['mask']
        sense = self.data[subject]['sense'] # [C, D, H, W]

        imspace = imspace[:, :, index] # [D, dyn, W]
        imspace_undersampled_mask = imspace_undersampled_mask[:, :, index]    
        sense = sense[:, :, index] # [C, D, W]
            
        imspace = imspace.transpose(1, 0, 2) # [T, D, W]
        imspace_undersampled_mask = imspace_undersampled_mask.transpose(1, 0, 2)
        
        imspace_undersampled_mask = imspace_undersampled_mask.astype(np.int8)
    
        sense = sense[:,:, np.newaxis] # [C, T, D, W]
        kspace = np.fft.fft(sense * imspace[np.newaxis], axis=-1)
        estimate = np.sum(np.fft.ifft(kspace, axis=-1) * sense.conj(), 0)
        data = {
            'target': np.ascontiguousarray(imspace), # [T, D, W]
            'sense': np.ascontiguousarray(sense), # [C, T, D, W]
            'mask': imspace_undersampled_mask[np.newaxis, ..., np.newaxis], # [1, T, D, W, 1]
            'measurements': kspace, # [C, T, D, W]
            'estimate': estimate # [T, D, W]
        }
        return data

    def augment_data(self, imspace, sense): # [T, D, W], [C, D, W]
        imspace, sense = self.random_crop(imspace, sense)
        if np.random.rand() < 0.5:
            imspace = np.flip(imspace, 2) 
            sense = np.flip(sense, 2)
        if np.random.rand() < 0.5:
            imspace = np.flip(imspace, 1) 
            sense = np.flip(sense, 1)
        if np.random.rand() < 0.5:
            imspace = np.flip(imspace, 0)
        imspace = np.roll(imspace, np.random.randint(0, imspace.shape[0]), 0)
        return imspace, sense

    def random_crop(self, imspace, sense):
        h, w = imspace.shape[-2], imspace.shape[-1]
        new_h, new_w = 25 if not h == 1 else 1, self.crop
        if h == new_h:
            top = 0
        else:
            top = np.random.randint(0, h - new_h)
        if w == new_w:
            left = 0
        else:
            left = np.random.randint(0, w - new_w)
        imspace = imspace[..., top:top + new_h, left:left + new_w]
        sense = sense[..., top:top + new_h, left:left + new_w]
        return imspace, sense


    







