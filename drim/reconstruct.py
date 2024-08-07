__author__ = 'Kai'

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

def reconstruct(config):

    train_config = configparser.ConfigParser()
    train_config.read(os.path.join(sys.argv[2], 'config.ini'))

    train_config = train_config['train']

    checkpoint = int(sys.argv[4])
       
    ndynamic = None
    # if config['ndynamic']:
    #     ndynamic = int(config['ndynamic'])
    
    dataset = MRData(
    sys.argv[3],
    int(config['crop']), 
    True)
    dataloader = data.DataLoader(
        dataset, int(config['batch-size']), drop_last=False, 
        sampler=data.SequentialSampler(dataset),
        num_workers=int(config['num-workers']))
    

    config['saved-model'] = os.path.join(
        sys.argv[2], 'network-parameters',
        'checkpoint' + sys.argv[4] + '.pt')
    network, initrim, gradrim = load_model(config, train_config, checkpoint)
    logger.info(f"Loaded network parameters from {config['saved-model']}.")


    if config.getboolean('randomize-bins'):
        logger.info("Using randomized bins...")
        with open(
            os.path.join(
                dataloader.dataset.data_path, 'imspace_header'
            ), 'rb'
        ) as f:
            nslice = pickle.load(f)[0][0]
        iterator = [
            None, 
            [torch.randperm(
                10, generator=torch.manual_seed(n)) for n in range(nslice)
            ]
        ]
    else:
        iterator = [None]

    for idx in iterator:
        with torch.no_grad():
            if not 'ndynamic' in config:
                reconstruct_data(
                    config, train_config, dataloader,
                    network, initrim, gradrim,
                    random_bin_idxs=idx)
            else:
                total_n_dynamic = dataset.data[
                    list(dataset.data.keys())[0]]['imspace'].shape[1]
                it_range = range(
                    0, total_n_dynamic, ndynamic)[:total_n_dynamic // ndynamic]
                for dynamic in it_range:
                    dataset.set_respiration(
                        dynamic, ndynamic, int(config['nbin']))
                    reconstruct_data(
                        config, train_config, dataloader,
                        network, initrim, gradrim,
                        random_bin_idxs=idx,
                        file_suffix=f'dynamics{dynamic}-{dynamic + ndynamic - 1}',
                        start_dynamic=dynamic
                    )


def reconstruct_data(
    config, train_config, dataloader, network, initrim, gradrim,
    random_bin_idxs=None, file_suffix='', start_dynamic=0
):
    start_time = time.time()
    
    if 'estimate' in config['wdir']:
        logger.info(
            f"Found the term 'estimate' in 'wdir' argument {config['wdir']}, "
            "meaning input will be saved and function be returned without "
            "reconstruction."
        )
        
        extracted_data = [
            dataloader.dataset[i] for i in range(len(dataloader.dataset))]
        
        estimate = np.stack([data['estimate'] for data in extracted_data], 2)
        with open(config['wdir'], 'wb') as f:
            pickle.dump(estimate, f, protocol=pickle.HIGHEST_PROTOCOL)
        if config.getboolean('save-target'):
            logger.info("Saving target, too.")
            target = np.stack([data['target'] for data in extracted_data], 2)
            with open(config['wdir'].replace('estimate', 'target'), 'wb') as f:
                pickle.dump(target, f, protocol=pickle.HIGHEST_PROTOCOL)
        return
        
         
    recons = []
    targets = []
    
    for batch in dataloader:

        target = torch.view_as_real(batch['target'].to(
            device=config['device'], dtype=torch.cfloat))
        estimate = torch.view_as_real(batch['estimate'].to(
            device=config['device'], dtype=torch.cfloat))
        measurements = torch.view_as_real(batch['measurements'].to(
            device=config['device'], dtype=torch.cfloat))
        sense = torch.view_as_real(batch['sense'].to(
            device=config['device'], dtype=torch.cfloat))
        mask = batch['mask'].to(device=config['device'], dtype=torch.bool)

        if not random_bin_idxs is None:
            for i, idx in enumerate(random_bin_idxs):
                newestimate = estimate[:, :, i].clone()
                estimate[:, idx, i] = newestimate
                newmeasurements = measurements[:, :, :, i].clone()
                measurements[:, :, idx, i] = newmeasurements
                newmask = mask[:, :, :, i].clone()
                mask[:, :, idx, i] = newmask

        if train_config.getboolean('autocast'):
            with autocast():
                hidden = initrim(estimate.moveaxis(-1, 1))
        else:
            hidden = initrim(estimate.moveaxis(-1, 1))
        for _ in range(int(config['niteration'])):
            gradient = gradrim(estimate, measurements, sense, mask)
            network_input = torch.cat(
                (estimate.moveaxis(-1, 1), gradient), 1) # [B, 4, T, D, W]
            if train_config.getboolean('autocast'):
                with autocast():
                    estimate_step, hidden = network(network_input, hidden)
            else:
                estimate_step, hidden = network(network_input, hidden)
            estimate = estimate + estimate_step.moveaxis(1, -1) # [B, T, D, W, 2]

        
        if not random_bin_idxs is None:
            for i, idx in enumerate(random_bin_idxs):
                inv_idx = torch.argsort(idx)
                newestimate = estimate[:, :, i].clone()
                estimate[:, inv_idx, i] = newestimate

        recons.append(estimate)
        if config.getboolean('save-target'):
            targets.append(target)
    if config.getboolean('save-target'):
        target = torch.cat(targets, 0).moveaxis(0, 3) # [T, D, H, W, 2]
    reconstruction = torch.cat(recons, 0).moveaxis(0, 3) # [T, D, H, W, 2]
 
    reconstruction = torch.view_as_complex(reconstruction).cpu().numpy()
    reconstruction = reconstruction.transpose((0, 3, 2, 1))
    print(reconstruction.shape)
    file_name = 'retroAItemp'
    mat_file = sys.argv[3] + file_name + '.mat'
    mat = loadmat(mat_file)
    print(np.array(mat['kData']).shape)
    print(reconstruction.shape)

    new_data = {'dData': np.abs(reconstruction)}
    mat.update(new_data)
    file_name = 'retroAItemp_DRIM'
    mat_file = sys.argv[3] + file_name + '.mat'
    savemat(mat_file, mat)
    end_time = time.time()
    print(end_time - start_time)
    exit()
    if reconstruction.shape[1] == 1:
        frames = []
        gif_name = 'reconstruction_model_10.gif'
        frame_duration = 30
        for timeframe in range(32):
            timeframe_slice = np.abs(reconstruction[timeframe, 0, :, :])
            image_space_normalized = cv2.normalize(timeframe_slice, None, 0, 255, cv2.NORM_MINMAX)
            frames.append(image_space_normalized) 
        imageio.mimsave(gif_name, frames, fps = frame_duration, loop=0)
    else:
        for slice in range(reconstruction.shape[1]):
            frames = []
            gif_name = 'gif-slices-checkerboard/drim_recon_slice_{}.gif'.format(slice)
            frame_duration = 10
            for timeframe in range(reconstruction.shape[0]):
                timeframe_slice = np.abs(reconstruction[timeframe, slice, :, :])
                image_space_normalized = cv2.normalize(timeframe_slice, None, 0, 255, cv2.NORM_MINMAX)
                frames.append(image_space_normalized)
            imageio.mimsave(gif_name, frames, fps = frame_duration, loop=0)
    
    logger.info("Finished reconstruction.")
    exit()

    if not config.getboolean('sorted'):
        dynslice = dataloader.dataset.dynamic_slice
        if config.getboolean('save-target'):
            tarslices = []
            iterator = zip(
                target.moveaxis(1, 0),
                reconstruction.moveaxis(1, 0),
                dataloader.dataset.data[''][
                    f'respiration-{dynslice.start}-{dynslice.stop - 1}'
                    f'dynamics-{dataloader.dataset.n_bin}bins'][0]
            )
        else:
            iterator = zip(
                reconstruction.shape[1] * [None],
                reconstruction.moveaxis(1, 0),
                dataloader.dataset.data[''][
                    f'respiration-{dynslice.start}-{dynslice.stop - 1}'
                    f'dynamics-{dataloader.dataset.n_bin}bins'][0]
            )
        recslices = []
        
        for tarslice, recslice, slicebreath in iterator: # [dyn, H, W, 2]
            if config.getboolean('save-target'):
                tardynamics = []
            recdynamics = []
            for bin in range(10):
                if not slicebreath[bin]:
                    if config.getboolean('save-target'):
                            tardynamics.append(torch.zeros_like(tarslice[0]))
                    recdynamics.append(torch.zeros_like(recslice[0]))
                else:
                    if config.getboolean('save-target'):
                        tardynamics.append(tarslice[
                            slicebreath[bin][
                                int((len(slicebreath[bin]) - 1)/2)
                            ][0]
                        ])
                    recdynamics.append(recslice[
                        slicebreath[bin][
                            int((len(slicebreath[bin]) - 1)/2)
                        ][0]
                    ])
            if config.getboolean('save-target'):
                tarslices.append(torch.stack(tardynamics)) # [T, H, W, 2]
            recslices.append(torch.stack(recdynamics))
        if config.getboolean('save-target'):
            target = torch.stack(tarslices, 1)
        reconstruction = torch.stack(recslices, 1) # [T, D, H, W, 2]
    
    if config.getboolean('save-target'):
        target = torch.view_as_complex(target).cpu().numpy()
    reconstruction = torch.view_as_complex(reconstruction).cpu().numpy()
    
    if config.getboolean('nifti'):
        logger.info(f"Saving reconstruction as nifti at {config['wdir']}.")
        reconstruction = nib.Nifti1Image(
            np.abs(reconstruction).real, np.eye(4))
        nib.save(reconstruction, config['wdir'])
        if config.getboolean('save-target'):
            logger.info("Saving target, too.")
            target = nib.Nifti1Image(np.abs(target).real, np.eye(4))
            nib.save(target, config['wdir'] + '_target')
    else:
        logger.info(f"Saving reconstruction as pickle at {config['wdir']}.")
        with open(config['wdir'] + file_suffix, 'wb') as f:
            pickle.dump(reconstruction, f, protocol=pickle.HIGHEST_PROTOCOL)
        if config.getboolean('save-target'):
            logger.info("Saving target, too.")
            with open(config['wdir'] + file_suffix + '_target', 'wb') as f:
                pickle.dump(target, f, protocol=pickle.HIGHEST_PROTOCOL)