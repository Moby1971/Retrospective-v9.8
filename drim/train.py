"""Train a model"""
__author__ = 'Kai'

import rim
from data_sampler import MRData
from validate import *

import os
import time
import math
import logging
import configparser
import matplotlib.pyplot as plt

import pandas as pd
import torch
from torch.utils.tensorboard import SummaryWriter
from torch.cuda.amp import GradScaler, autocast
import torch.utils.data as data

import logging
logger = logging.getLogger(__name__)

def time_since(since):
    """Get the time elapsed since since"""
    now = time.time()
    s = now - since
    m = math.floor(s / 60)
    t = math.floor(m / 60)
    s -= m * 60
    m -= t * 60
    return t, m, s

def load_model(config):
    old_config = configparser.ConfigParser()
    old_config.read(os.path.join(
        os.path.dirname(os.path.dirname(config['saved-model'])), 'config.ini')
    )
    old_config = old_config['train']
    nfeature = int(old_config['nfeature'])
    kernel = []
    for kern in old_config['kernel'].split():
        if kern == 'None':
            kernel.append(None)
        else:
            kernel.append([int(k) for k in kern])
    network = rim.RecurrentInferenceMachine(
        nfeature=nfeature, kernel=kernel,
        temporal_rnn=config.getboolean('temporal-rnn'), mute=False)
    initrim = rim.InitRim(2, 2 * [nfeature], kernel, mute=False)
    
    if 'fourier-dim' in old_config:
        fourier_dim = [int(d) for d in old_config['fourier-dim'].split()]
    else:
        fourier_dim = -1
    gradrim = rim.GradRim(fourier_dim=fourier_dim)

    load = torch.load(config['saved-model'], 
        map_location=lambda storage, loc: storage.cpu())
        
    network.load_state_dict(load['rim'])
    try:
        initrim.load_state_dict(load['initrim'])
    except:
        initrim.load_state_dict(load['init'])
    network = network.to(device=config['device'])
    initrim = initrim.to(device=config['device'])
    
    optimizer = torch.optim.Adam(
        list(network.parameters()) + list(initrim.parameters()), 
        lr=float(old_config['lr']))
    optimizer.load_state_dict(load['optimizer'])
    
    step = int(config['saved-model'].split(
        'checkpoint')[1].split('.pt')[0]) + 1
    
    return network, initrim, gradrim, optimizer, step

def train_model(config):
    logger = logging.getLogger(__name__)

    if 'saved-model' in config:
        network, initrim, gradrim, optimizer, step = load_model(config)
    else:
        nfeature = int(config['nfeature'])
        kernel = []
        for kern in config['kernel'].split():
            if kern == 'None':
                kernel.append(None)
            else:
                kernel.append([int(k) for k in kern])
        network = rim.RecurrentInferenceMachine(
            nfeature=nfeature, kernel=kernel,
            temporal_rnn=config.getboolean('temporal-rnn'), mute=False)
        initrim = rim.InitRim(2, 2 * [nfeature], kernel, mute=False)
        network = network.to(device=config['device'])
        initrim = initrim.to(device=config['device'])
        optimizer = torch.optim.Adam(
            list(network.parameters()) + list(initrim.parameters()), 
            lr=float(config['lr']))
        gradrim = rim.GradRim(
            fourier_dim=[int(d) for d in config['fourier-dim'].split()])
        step = 0

    network.train()
    initrim.train()
    gradrim.train()
    compute_loss = torch.nn.L1Loss(reduction='mean')

    scheduler = torch.optim.lr_scheduler.MultiStepLR( 
        optimizer, 
        [int(ms) for ms in config['milestones'].split()], 
        gamma=float(config['gamma'])
    )

    dataset = MRData(
        config['data-dir'],
        config['undersampled_data_dir'],
        int(config['crop']), 
        True)
    dataloader = data.DataLoader(
        dataset, int(config['batch-size']), shuffle=True,
        drop_last=True, num_workers=int(config['num-workers']))
    
    dataloader_visualize = data.DataLoader(
        dataset, int(config['batch-size']), shuffle=False,
        drop_last=True, num_workers=int(config['num-workers']))
    
    loss_weights = torch.logspace(
        -1, 0, steps=int(config['niteration'])).to(next(network.parameters()))

    logdir = os.path.join(config['train-dir'], 'logs', 'tensorboard')
    writer = SummaryWriter(log_dir=logdir)

    traindir = os.path.join(config['train-dir'], 'network-parameters')

    mixed_precision_scaler = GradScaler()

    #VISUALIZE ESTIMATE
    viz = []
    for batch in dataloader_visualize:
        viz.append(batch['estimate'])
    viz = torch.cat(viz, 0).moveaxis(0, 3)
    print(viz.shape)
    plot_images(config, '', viz, 'train', 'undersampled', 0, writer )
    #viz = torch.view_as_complex(viz).cpu().numpy()
    # plt.imshow(np.abs(viz[0, 0, :, :]), cmap='gray')
    # plt.show()


    running_loss = .0
    start = time.time()
    logger.info(f'Starting training...')
    for epoch in range(int(step / len(dataloader)), int(config['nepoch'])):
        recons = []
        for batch in dataloader:
            
            

            optimizer.zero_grad()            
            target = torch.view_as_real(batch['target'].to(
                device=config['device'], dtype=torch.cfloat))
            estimate = torch.view_as_real(batch['estimate'].to(
                device=config['device'], dtype=torch.cfloat))
            measurements = torch.view_as_real(batch['measurements'].to(
                device=config['device'], dtype=torch.cfloat))
            sense = torch.view_as_real(batch['sense'].to(
                device=config['device'], dtype=torch.cfloat))
            mask = batch['mask'].to(device=config['device'], dtype=torch.bool)

            
               

            loss = torch.Tensor([0.]).to(
                device=config['device'], dtype=torch.float)
            
            if config.getboolean('autocast'):
                with autocast():
                    hidden = initrim(estimate.moveaxis(-1, 1))
            else:
                hidden = initrim(estimate.moveaxis(-1, 1))

            for iteration, weight in zip(
                range(int(config['niteration'])), loss_weights):
                gradient = gradrim(estimate, measurements, sense, mask) # [B, T, D, W]
                network_input = torch.cat(
                    (estimate.moveaxis(-1, 1), gradient), 1) # [B, 4, T, D, W]
                if config.getboolean('autocast'):
                    with autocast():
                        estimate_step, hidden = network(network_input, hidden) # [B, 2, T, D, W]
                else:
                    estimate_step, hidden = network(network_input, hidden)
                estimate = estimate + estimate_step.moveaxis(1, -1)

                
                it_loss = compute_loss(estimate, target)
                loss = loss + weight * it_loss
                if iteration % int(config['truncate']) == int(
                    config['truncate']) - 1:
                    loss = loss / int(config['truncate'])
                    mixed_precision_scaler.scale(loss).backward()
                    loss = torch.Tensor([0.]).to(
                        device=config['device'], dtype=torch.float)
                    hidden = [h.detach() for h in hidden]
                    estimate = estimate.detach()
            recons.append(estimate)
            
            mixed_precision_scaler.step(optimizer)
            mixed_precision_scaler.update()
            running_loss += it_loss.item()
            writer.add_scalar(
                'Train/Loss/FinalLoss', it_loss.item(), step + 1)
            
            if step % int(config['print-freq']
            ) == int(config['print-freq']) - 1:
                logger.info('[{:d}, {:7,d}] train loss: {:e}.'.format(
                    epoch + 1,
                    step + 1,
                    running_loss / int(config['print-freq'])
                ))
                running_loss = .0
            if step % int(config['checkpoint-freq']
            ) == int(config['checkpoint-freq']) - 1:
                torch.save(
                    {
                        'rim': network.state_dict(), 
                        'initrim': initrim.state_dict(),
                        'optimizer': optimizer.state_dict()
                    },
                    os.path.join(
                        traindir, f'checkpoint{str(step + 1)}.pt'
                    )
                )
                timesincestart = time_since(start)
                logger.info(
                    'Checkpoint made at epoch {}, iteration {} | '
                    'Time elapsed since start: {}'.format(
                        epoch + 1, step + 1, '%dt %dm %ds' % timesincestart))
                if timesincestart[
                    0] + timesincestart[1]/60. > float(config['time-limit']):
                    print('Time limit exceeded, aborting train function...')
                    return
            step += 1
        scheduler.step()
        #VISUALIZE AFTER EPOCH TRAINING ESTIMATES
        # reconstruction = torch.cat(recons, 0).moveaxis(0, 3)
        # reconstruction = torch.view_as_complex(reconstruction).cpu().numpy()
        # plt.imshow(np.abs(reconstruction[0, 0, :, :]), cmap='gray')
        # plt.show()

        #VISUALIZE AFTER EPOCH LATEST CHECKPOINT

        torch.cuda.empty_cache()
        network.eval()
        initrim.eval()
        gradrim.eval()
        with torch.no_grad():
        
            recons_visualize = []
            for batch in dataloader_visualize:
                        
                target = torch.view_as_real(batch['target'].to(
                    device=config['device'], dtype=torch.cfloat))
                estimate = torch.view_as_real(batch['estimate'].to(
                    device=config['device'], dtype=torch.cfloat))
                measurements = torch.view_as_real(batch['measurements'].to(
                    device=config['device'], dtype=torch.cfloat))
                sense = torch.view_as_real(batch['sense'].to(
                    device=config['device'], dtype=torch.cfloat))
                mask = batch['mask'].to(device=config['device'], dtype=torch.bool)

                
                if config.getboolean('autocast'):
                    with autocast():
                        hidden = initrim(estimate.moveaxis(-1, 1))
                else:
                    hidden = initrim(estimate.moveaxis(-1, 1))

                for iteration, weight in zip(
                    range(int(config['niteration'])), loss_weights):
                    gradient = gradrim(estimate, measurements, sense, mask) # [B, T, D, W]
                    network_input = torch.cat(
                        (estimate.moveaxis(-1, 1), gradient), 1) # [B, 4, T, D, W]
                    if config.getboolean('autocast'):
                        with autocast():
                            estimate_step, hidden = network(network_input, hidden) # [B, 2, T, D, W]
                    else:
                        estimate_step, hidden = network(network_input, hidden)
                    estimate = estimate + estimate_step.moveaxis(1, -1)

                    
                recons_visualize.append(estimate)
            reconstruction_viz = torch.cat(recons_visualize, 0).moveaxis(0,3)
            reconstruction_viz = torch.view_as_complex(reconstruction_viz)
            plot_images(config, '', reconstruction_viz, 'train', 'undersampled', epoch, writer)