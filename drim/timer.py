import timeit
import torch
import numpy as np
from torch.cuda.amp import autocast
import logging
import rim

def time_model(config):
    logger = logging.getLogger(__name__)
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
    initrim = rim.InitRim(
        2, 2 * [nfeature], kernel, mute=False)
    network = network.to(device=config['device'])
    initrim = initrim.to(device=config['device'])
    gradrim = rim.GradRim(
        fourier_dim=[int(d) for d in config['fourier-dim'].split()])
    network.eval()

    ncoil = int(config['ncoil'])
    ndynamic = int(config['ndynamic'])
    nlocation = int(config['nlocation'])
    width = int(config['width'])
    batch = int(config['height'])
    
    mask = torch.tensor(np.random.choice([True, False],
        size=(batch, 1, ndynamic, nlocation, width, 1))).cuda()
    estimate = torch.rand(
        batch, ndynamic, nlocation, width, 2, dtype=torch.float).cuda()
    measurements = torch.rand(
        batch, ncoil, ndynamic, nlocation, width, 2, dtype=torch.float).cuda()
    sense = torch.rand(
        batch, ncoil, 1, nlocation, width, 2, dtype=torch.float).cuda()
    gradient = torch.zeros(batch, 2, ndynamic, nlocation, width).to(estimate)

    def inference():
        est = estimate
        grad = gradient
        with autocast():
            hidden = initrim(est.moveaxis(-1, 1))
        for iteration in range(int(config['niteration'])):
            with autocast():
                if iteration != 0:
                    grad = gradrim(est, measurements, sense, mask)
                network_input = torch.cat(
                    (est.moveaxis(-1, 1), grad), 1) # [B, 4, T, D, W]
                estimate_step, hidden = network(network_input, hidden)
                est = est + estimate_step.moveaxis(1, -1)

    with torch.no_grad():
        logger.info(
            timeit.repeat(
                'inference()',
                globals=locals(),
                number=int(config['nbatch']),
                repeat=int(config['nrepeat']))
        )