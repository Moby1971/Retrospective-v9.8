__author__ = 'BMEP'

import sys
sys.path.insert(0, './drim')
import os
import configparser
import pickle

import numpy as np
import nibabel as nib

import torch
import torch.utils.data as data
from torch.cuda.amp import autocast
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
import imageio
import cv2

from data_sampler import MRData
from validate import load_model
import skimage.metrics as measure
import numpy.fft as fft

import logging
import time
logger = logging.getLogger(__name__)

#print(torch.backends.mps.is_available())
#print(torch.backends.mps.is_built())

# sys.argv[1] = reconstruct
# sys.argv[2] = train directory
# sys.argv[3] = temporary directory where mat files are located
# sys.argv[4] = checkpoint 
# sys.argv[5] = drim config file

# print(sys.argv[1])
# print(sys.argv[2])
# print(sys.argv[3])
# print(sys.argv[4])
# print(sys.argv[5])


def reconstruct(config):

    train_config = configparser.ConfigParser()
    train_config.read(os.path.join(sys.argv[2], 'config.ini'))
    train_config = train_config['train']
    checkpoint = int(sys.argv[4])
  
    dataset = MRData(sys.argv[3])

    dataloader = data.DataLoader(
        dataset, int(config['batch-size']), drop_last=False, 
        sampler=data.SequentialSampler(dataset),
        num_workers=int(config['num-workers']))
  
    config['saved-model'] = os.path.join(
        sys.argv[2], 'network-parameters',
        'checkpoint' + sys.argv[4] + '.pt')

    network, initrim, gradrim = load_model(config, train_config, checkpoint)
    logger.info(f"Loaded network parameters from {config['saved-model']}.")

    with torch.no_grad():
        reconstruct_data(config, train_config, dataloader, network, initrim, gradrim, sys.argv[3])




def reconstruct_data(config, train_config, dataloader, network, initrim, gradrim, temp_dir):
    start_time = time.time()

  # for param in network.parameters():
  #      print(param)
  #      break
         
    recons = []
    
    for batch in dataloader:

        estimate = torch.view_as_real(batch['estimate'].to(device=config['device'], dtype=torch.cfloat))
        measurements = torch.view_as_real(batch['measurements'].to(device=config['device'], dtype=torch.cfloat))
        sense = torch.view_as_real(batch['sense'].to(device=config['device'], dtype=torch.cfloat))
        mask = batch['mask'].to(device=config['device'], dtype=torch.bool)
       
        #print(f'sense shape: {sense.shape}')
        #print(f'mask shape: {mask.shape}')
        #print(f'measurement shape: {measurements.shape}')
        #print(f'estimate shape: {estimate.shape}')

        #print(estimate[0,0,0,:,0])

        hidden = initrim(estimate.moveaxis(-1, 1))
        for _ in range(int(config['niteration'])):
            gradient = gradrim(estimate, measurements, sense, mask)
            network_input = torch.cat((estimate.moveaxis(-1, 1), gradient), 1) # [B, 4, T, D, W]
            estimate_step, hidden = network(network_input, hidden)
            estimate = estimate + estimate_step.moveaxis(1, -1) # [B, T, D, W, 2]
     
        recons.append(estimate)


    reconstruction = torch.cat(recons, 0).moveaxis(0, 3) # [T, D, H, W, 2]
    reconstruction = torch.view_as_complex(reconstruction).cpu().numpy()
    reconstruction = reconstruction.transpose((0, 3, 2, 1))
    
    file_name = 'retroAItemp'
    mat_file = sys.argv[3] + file_name + '.mat'
    mat = loadmat(mat_file)

    new_data = {'aiData': np.abs(reconstruction)}
    mat.update(new_data)
    print(reconstruction.shape)
    print(np.array(mat['aiData']).shape)

    file_name = 'retroAItemp_DRIM'
    mat_file = sys.argv[3] + file_name + '.mat'
    savemat(mat_file, mat)

    end_time = time.time()
    print(end_time - start_time)
   
    logger.info("Finished reconstruction.")
    exit()
