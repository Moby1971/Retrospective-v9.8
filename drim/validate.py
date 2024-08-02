__author__ = 'Kai'

import sys
sys.path.insert(0, './drim')
import os
import configparser
import pickle

import torch
import torch.utils.data as data
from torch.cuda.amp import autocast

from torch.utils.tensorboard import SummaryWriter
import torchvision

import torchmetrics

import numpy as np
import pandas as pd
import skimage.metrics as measure
from interval import Interval

import rim
from data_sampler import MRData

import logging
logger = logging.getLogger(__name__)

def load_model(config, train_config, checkpoint):
    nfeature = int(train_config['nfeature'])
    kernel = []
    for kern in train_config['kernel'].split():
        if kern == 'None':
            kernel.append(None)
        else:
            kernel.append([int(k) for k in kern])
    network = rim.RecurrentInferenceMachine(
        nfeature=nfeature, kernel=kernel,
        temporal_rnn=train_config.getboolean('temporal-rnn'))
    initrim = rim.InitRim(2, 2 * [nfeature], kernel)
    
    if 'fourier-dim' in train_config:
        fourier_dim = [int(d) for d in train_config['fourier-dim'].split()]
    else:
        fourier_dim = -1
    gradrim = rim.GradRim(fourier_dim=fourier_dim)

    load = torch.load(
        os.path.join(
            sys.argv[2],
            'network-parameters',
            f"checkpoint{checkpoint}.pt"
        ),
        map_location=lambda storage, loc: storage.cpu()
    )
    network.load_state_dict(load['rim'])
    initrim.load_state_dict(load['initrim'])

    network = network.to(device=config['device'])
    initrim = initrim.to(device=config['device'])
    
    return network, initrim, gradrim

def validate_model(config):

    train_config = configparser.ConfigParser()
    train_config.read(os.path.join(sys.argv[2], 'config.ini'))
    train_config = train_config['train']

    #Default to validating the whole range of checkpoints
    if not config['min-max-checkpoint']:
        checkpoints = sorted(filter(lambda x: 'checkpoint' in x,
            os.listdir(os.path.join(
                sys.argv[2], 'network-parameters'))), 
            key=lambda x: int(x.split('checkpoint')[1].split('.')[0]))
        minmaxcheckpoint = (
            int(checkpoints[0].split('checkpoint')[1].split('.')[0]), 
            int(checkpoints[-1].split('checkpoint')[1].split('.')[0]) + int(
                train_config['checkpoint-freq']))
    else:
        minmaxcheckpoint = [
            int(mmc) for mmc in config['min-max-checkpoint'].split()]

    
    #Write columns of csv log file if none exists yet
    if config['test']:
        assert config['ndynamic']
        ndynamic = int(config['ndynamic'])
        if os.path.exists(os.path.join(
            sys.argv[2], 'logs', f"{config['test']}.csv")):
            tmp_edition = '1'
            while True:
                tmp_edition = str(int(tmp_edition) + 1)
                csvname = os.path.join(
                    sys.argv[2], 'logs',
                    config['test'] + f'_({tmp_edition}).csv')
                if not os.path.exists(csvname):
                    config['test'] = config['test'] + f'_({tmp_edition})'
                    break
        else:
            csvname = os.path.join(
                sys.argv[2], 'logs', f"{config['test']}.csv")
    else:
        csvname = os.path.join(sys.argv[2], 'logs', 'validation.csv')
        ndynamic = None
    if not os.path.exists(csvname):
        pd.DataFrame(
            columns=[
                'Checkpoint', 'T', 'Mask', 'Metric', 'Score',
                'Data', 'Number of gaps', 'Temporal gaps',
                'Dynamics', '4d sorted', 'Number of bins']
            ).to_csv(csvname, header=True, index=False, mode='w')
    writer = SummaryWriter(
        log_dir=os.path.join(sys.argv[2], 'logs', 'tensorboard'))
    logger.info(f'Writing results to {csvname}')
    nbin = int(config['nbin'])

    dataset = MRData(
    sys.argv[3],
    config['undersampled_data_dir'],
    int(config['crop']), 
    True)
    dataloader = data.DataLoader(
        dataset, 
        int(config['batch-size']), 
        drop_last=False, 
        sampler=data.SequentialSampler(dataset),
        num_workers=int(config['num-workers']))
    
    subject = os.path.basename(sys.argv[3])
    
    for checkpoint in map(str, range(*minmaxcheckpoint, int(config['checkpoint-freq']) * int(train_config['checkpoint-freq']))):
        print(checkpoint)
        network, initrim, gradrim = load_model(config, train_config, checkpoint)

        if config.getboolean('randomize-bins'):
            with open(
                os.path.join(
                    dataloader.dataset.data_path, 'imspace_header'
                ), 'rb'
            ) as f:
                nslice = pickle.load(f)[0][0]
            iterator = [
                None, 
                [torch.randperm(
                    nbin, generator=torch.manual_seed(n)) for n in range(nslice)
                ]
            ]
        else:
            iterator = [None]

        for idx in iterator:
            with torch.no_grad():
                if not config['ndynamic']:
                    validate_dataset(
                        config, train_config, subject, dataloader,
                        checkpoint, network, initrim, gradrim, csvname, writer,
                        random_bin_idxs=idx)
                else:
                    # Loop through whole acquisition sequence in batches of 
                    # FLAGS.n_dynamic length, except for last batch if it is of
                    # shorter length.
                    total_n_dynamic = dataset.data[
                        list(dataset.data.keys())[0]]['imspace'].shape[1]
                    it_range = range(
                        0, total_n_dynamic, int(config['ndynamic'])
                    )[:total_n_dynamic // int(config['ndynamic'])]
                    for dynamic in it_range:
                        dataloader.dataset.set_respiration(
                            dynamic,
                            int(config['ndynamic']),
                            int(config['nbin'])
                        )
                        validate_dataset(
                            config, train_config, subject,
                            dataloader, checkpoint, network,
                            initrim, gradrim, csvname, writer,
                            random_bin_idxs=idx
                        )

def validate_dataset(
    config, train_config, subject, dataloader, checkpoint,
    network, initrim, gradrim, csvname, writer, random_bin_idxs=None):
    
    metrics = {
        'nrmse': nrmse,
        'ssim2d': lambda y, x: SSIM()(y, x).item(),
        'ssim4d': lambda y, x: SSIM(dim=4)(y, x).item(),
        'ms-ssim': lambda y, x: ms_ssim(y, x).item(),
        'ssim': lambda y, x: ssim(y, x, win_size=int(config['ssim-win-size'])),
        'psnr': psnr,
        'tenengrad-target': tenengrad,
        'tenengrad-reconstruction': lambda y, x: tenengrad(x, y),
        'gradient-entropy-target': gradient_entropy,
        'gradient-entropy-reconstruction': lambda y, x: gradient_entropy(x, y)
    }

    recons = []
    targets = []
    viz = []
    trgts = []
    for batch in dataloader:

        target = torch.view_as_real(batch['target'].to(
            device=config['device'], dtype=torch.cfloat)) # [B/H, T, D, W, 2]
        estimate = torch.view_as_real(batch['estimate'].to(
            device=config['device'], dtype=torch.cfloat))
        measurements = torch.view_as_real(batch['measurements'].to(
            device=config['device'], dtype=torch.cfloat))
        sense = torch.view_as_real(batch['sense'].to(
            device=config['device'], dtype=torch.cfloat))
        mask = batch['mask'].to(device=config['device'], dtype=torch.bool)

        
        viz.append(batch['estimate'])
        trgts.append(batch['target'])
        

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
        for iteration in range(int(config['niteration'])):
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
        targets.append(target)

    viz = torch.cat(viz, 0).moveaxis(0, 3)
    plot_images(config, '', viz, 'Original data', 'reconstruction', 0, writer )
    trgts = torch.cat(trgts, 0).moveaxis(0, 3)
    plot_images(config, '', trgts, 'Original data', 'target', 0, writer, 'target-gif' )

    target = torch.cat(targets, 0).moveaxis(0, 3) # [T, D, H, W, 2]
    reconstruction = torch.cat(recons, 0).moveaxis(0, 3) # [T, D, H, W, 2]

  

    if not random_bin_idxs is None:
        subject = f'{subject}-random-bin'
    if not hasattr(dataloader.dataset, 'dynamic_slice'):
        dyns = 120
    else:
        dynrange = dataloader.dataset.dynamic_slice
        dynrange = f'{dynrange.start}-{dynrange.stop - 1}'
        subject = f'{subject}-dynamics{dynrange}'
        dyns = int(config['ndynamic'])

    plot_images(
        config, torch.view_as_complex(target),
        torch.view_as_complex(reconstruction), subject,
        'reconstruction', checkpoint, writer) #dataloader.dataset.maskname
    
    if config['test']:
        gap_mask = torch.zeros(reconstruction.shape[:2], dtype=bool)
        gap_neighbor_mask = torch.zeros(reconstruction.shape[:2], dtype=bool)
        rec_gaps = [] # Only gaps, unaltered reconstructions
        tar_gaps = [] # Only gaps, selected from best fitting slice from all dynamics
        rec_neighbors = [] # Only unaltered reconstructions of gap neibors
        tar_neighbors = [] # Only gap neighbors
        rec_gap_removed = [] # Only non-gaps, stacked in slice-bin order
        tar_gap_removed = [] # Only non-gaps, stacked in slice-bin order
        rec_gap_removed_2 = [] # Only non-gaps, stacked in bin-slice order
        tar_gap_removed_2 = [] # Only non-gaps, stacked in bin-slice order
        rec_gaps_neighbors_removed = []
        tar_gaps_neighbors_removed = []
        rec_interpolated_gaps = [] # Only gaps, interpolated from neighbor reconstructions
        rec_interpolated = torch.zeros_like(
            reconstruction).copy_(reconstruction) # Full volume with interpolations
        tar_filled_gaps = torch.zeros_like(target).copy_(target) # Full volume replaced where possible, otherwise interpolated

        data = dataloader.dataset.data['']
        breathing_data, bin_edges, neighbor_idx = data[
            f'respiration-{dynrange}dynamics-{int(config["nbin"])}bins'
        ]
        breathing_data_total, _, _ = data[
            f'respiration-{int(config["nbin"])}bins'
        ]
        
        # Make masks locating temporal gaps and their neighboring slices
        gap_idx = [[], []]
        for slice_position, breathing_bins in enumerate(breathing_data):
            for bin in range(int(config['nbin'])):
                if breathing_bins[bin]:
                    rec_gap_removed.append(reconstruction[bin, slice_position])
                    tar_gap_removed.append(target[bin, slice_position])
                else:
                    gap_idx[0].append(bin)
                    gap_idx[1].append(slice_position)
                    gap_mask[bin, slice_position] = True
                    gap_neighbor_mask[bin, slice_position] = True
                    ubin = (bin - 1) % int(config['nbin'])
                    dbin = (bin + 1) % int(config['nbin'])
                    gap_neighbor_mask[dbin, slice_position] = True
                    gap_neighbor_mask[ubin, slice_position] = True
                    if slice_position:
                        gap_neighbor_mask[bin, slice_position - 1] = True
                    if  slice_position != reconstruction.shape[1] - 1:
                        gap_neighbor_mask[bin, slice_position + 1] = True
        for bin in range(int(config['nbin'])):
            for slice_position, breathing_bins in enumerate(breathing_data):
                if breathing_bins[bin]:
                    rec_gap_removed_2.append(
                        reconstruction[bin, slice_position])
                    tar_gap_removed_2.append(target[bin, slice_position])

        bin_dict = {}
        for bin_nr in range(1, int(config['nbin']) // 2):
            bin_dict[bin_nr] = lambda x: (x[bin_nr], x[bin_nr + 1])
            bin_dict[int(config['nbin']) - bin_nr
                ] = lambda x: (x[bin_nr], x[bin_nr + 1])
        bin_dict[0] = lambda x: (x[0], x[1])
        bin_dict[int(config['nbin']) // 2] = lambda x: (
            x[int(config['nbin']) // 2], x[int(config['nbin']) // 2 + 1])

        def replace_nearest_bin_target(bin, slice_position):
            edges = bin_dict[bin](bin_edges)
            bin_center = sum(edges) / 2.
            edges = Interval(edges[0], edges[1])
            candidates = [[], []]
            for b in ( 
                (bin - 2) % int(config['nbin']), 
                (bin - 1) % int(config['nbin']),
                bin,
                (bin + 1) % int(config['nbin']),
                (bin + 2) % int(config['nbin'])
            ):
                if breathing_data_total[slice_position][b]:
                    for dynamic, signal in breathing_data_total[
                        slice_position][b]:
                        if signal in edges:
                            candidates[0].append(dynamic)
                            candidates[1].append(
                                np.abs(signal - bin_center))
            if candidates[0]:
                closest = candidates[0][np.argmin(candidates[1])]
                replacement_slice = torch.view_as_real(
                    torch.tensor(
                        dataloader.dataset.data[
                        '']['imspace'][slice_position, closest]
                    ).moveaxis(1, 0)
                )
                return torch.tensor(replacement_slice).to(target)

        def interpolate(reconstruction, neighbors):
            interpol = torch.zeros_like(reconstruction[0, 0])
            if neighbors:
                for bin, slice_position in neighbors:
                    interpol += reconstruction[bin, slice_position]
                return interpol / len(neighbors)

        
        for s in range(len(breathing_data)):
            for b in range(int(config['nbin'])):
                if not gap_neighbor_mask[b, s]:
                    rec_gaps_neighbors_removed.append(reconstruction[b, s])
                    tar_gaps_neighbors_removed.append(target[b, s])
                elif not gap_mask[b, s]:
                    rec_neighbors.append(reconstruction[b, s])
                    tar_neighbors.append(target[b, s])
                else:
                    interpolation = interpolate(
                        reconstruction, neighbor_idx[(b, s)])
                    if interpolation is not None:
                        rec_interpolated[b, s] = interpolation
                    replacement_slice = replace_nearest_bin_target(b, s)
                    if replacement_slice is not None:
                        rec_gaps.append(reconstruction[b, s])
                        tar_gaps.append(replacement_slice)
                        rec_interpolated_gaps.append(rec_interpolated[b, s])
                        tar_filled_gaps[b, s] = replacement_slice
                    else:
                        tar_filled_gaps[b, s] = torch.zeros_like(target[0, 0])
                        rec_interpolated[b, s] = torch.zeros_like(
                            reconstruction[0, 0])
        
        n_gap = len(gap_idx[0])
        rec_gap_removed = torch.stack(rec_gap_removed).unsqueeze(0)
        tar_gap_removed = torch.stack(tar_gap_removed).unsqueeze(0)
        rec_gap_removed_2 = torch.stack(rec_gap_removed_2).unsqueeze(0)
        tar_gap_removed_2 = torch.stack(tar_gap_removed_2).unsqueeze(0)
        rec_gaps_neighbors_removed = torch.stack(
            rec_gaps_neighbors_removed).unsqueeze(0)
        tar_gaps_neighbors_removed = torch.stack(
            tar_gaps_neighbors_removed).unsqueeze(0)
        to_evaluate = []
        if n_gap > 0:
            """if rec_gaps:
                rec_gaps = torch.stack(rec_gaps).unsqueeze(0)
                tar_gaps = torch.stack(tar_gaps).unsqueeze(0)
                rec_neighbors = torch.stack(rec_neighbors).unsqueeze(0)
                tar_neighbors = torch.stack(tar_neighbors).unsqueeze(0)
                to_evaluate.extend((
                    ('gaps only - left as output', rec_gaps, tar_gaps),
                    ))"""
            if rec_interpolated_gaps and (
                isinstance(tar_gaps, list) and tar_gaps):
                rec_interpolated_gaps = torch.stack(
                    rec_interpolated_gaps).unsqueeze(0)
                to_evaluate.extend((
                    ('whole volume - gaps interpolated',
                     rec_interpolated, tar_filled_gaps),))
                    #('gaps only - interpolated', 
                    # rec_interpolated_gaps, tar_gaps),
                    #('gap neighbors only', rec_neighbors, tar_neighbors)))
        to_evaluate.extend([
            (
                'whole volume - gaps left as output', 
                reconstruction, tar_filled_gaps
            ),
            (
                'gaps removed - stacked in slice-bin order',
                rec_gap_removed, tar_gap_removed
            ),
            (
                'gaps removed - stacked in bin-slice order',
                rec_gap_removed_2, tar_gap_removed_2
            ),
            #(
            #    'gaps and gap neighbors removed',
            #    rec_gaps_neighbors_removed, tar_gaps_neighbors_removed
            #)
        ])
    else:
        to_evaluate = [('N/A', reconstruction, target)]
        n_gap = 0

    is_sorted = "4d-sorted" if config.getboolean("sorted") else "unsorted"
    for metric in config['metrics'].split():
        for label, r, t in to_evaluate:
            if config.getboolean('perslicescore'):
                for count, (r_, t_) in enumerate(zip(
                    r.reshape(-1, *r.shape[-2:]), t.reshape(-1, *t.shape[-2:]))
                ):
                    score = metrics[metric](t_, r_)
                    df = pd.DataFrame(
                        data={
                            'T': int(config['niteration']),
                            'Mask': dataloader.dataset.maskname, 
                            'Metric': metric, 
                            'Score': score, 
                            'Data': subject,
                            'Number of gaps': n_gap,
                            'Temporal gaps': label,
                            'Dynamics': dyns,
                            '4d sorted': config.getboolean('sorted'),
                            'Number of bins': int(config['nbin']),
                            'Slice': count % t.shape[1],
                            'Bin': count // t.shape[1],
                            },
                        index=[int(checkpoint)])
                    df.to_csv(csvname, header=False, mode='a')
            else:
                score = metrics[metric](t, r)
                logger.info(
                    f"Checkpoint {checkpoint} {subject} "
                    f"mask {metric} with {n_gap} " #{dataloader.dataset.maskname}
                    f"gaps and {int(config['nbin'])} bins; {label}: {score}")
                df = pd.DataFrame(
                    data={
                        'T': int(config['niteration']),
                        'Mask': 'mask', #dataloader.dataset.maskname
                        'Metric': metric, 
                        'Score': score, 
                        'Data': subject,
                        'Number of gaps': n_gap,
                        'Temporal gaps': label,
                        'Dynamics': dyns,
                        '4d sorted': config.getboolean('sorted'),
                        'Number of bins': int(config['nbin'])},
                    index=[int(checkpoint)]
                )
                df.to_csv(csvname, header=False, mode='a')
                writer.add_scalar(
                    f'Validation/{subject}/{is_sorted}/'
                    f'mask/{n_gap}_' #{dataloader.dataset.maskname}
                    f'{"_".join(label.split())}/{int(config["nbin"])}/{metric}',
                    score,
                    int(checkpoint)
                )

def nrmse(y, x): # [T, D, H, W]
    output = torch.sqrt(torch.sum(
        torch.pow(y - x, 2.)) / torch.sum(torch.pow(y, 2.)))
    return output.item()
    
def ssim(y, x, win_size=7): # [T, D, H, W]
    target = y.squeeze().cpu().numpy()
    recon = x.squeeze().cpu().numpy() # [T, D, H, W, 2]
    return measure.structural_similarity(
        target, recon, channel_axis=-1,
        data_range=target.max(), win_size=win_size)

def psnr(y, x): # [T, D, H, W]
    target = torch.abs(torch.view_as_complex(y)).squeeze().cpu().numpy()
    recon = torch.abs(torch.view_as_complex(x)).squeeze().cpu().numpy() # [T, D, H, W, 2]
    return measure.peak_signal_noise_ratio(
        target, recon, data_range=target.max())

def gradient_entropy(y, x): # [T, D, H, W, 2] or [D, H, W, 2]
    image = torch.abs(torch.view_as_complex(
        y)).squeeze().cpu().numpy() # [T, D, H, W] or [D, H, W]
    gradients = np.gradient(image.reshape(-1, image.shape[-1]), axis=1) # [TDH, W]
    abs_grad = np.abs(gradients)
    sum_abs_grad = np.sum(abs_grad, axis=1, keepdims=True) # [TDH, 1]
    norm_grad = abs_grad / sum_abs_grad
    norm_grad[np.where(norm_grad == np.nan)] = 0
    log = np.log2(norm_grad)
    log[np.where(log == -np.inf)] = 0
    entropy = - np.sum(norm_grad * log, axis=1) # [TDH]
    return entropy.mean()

def tenengrad(y, x): # [T, D, H, W, 2] or [D, H, W, 2]
    image = torch.abs(torch.view_as_complex(
        y)).squeeze().cpu().numpy() # [T, D, H, W] or [D, H, W]
    gradients = np.gradient(image.reshape(-1, *image.shape[-2:]), axis=(1, 2)) # 2 * [TD, H, W]
    tenengrad = np.sum(np.abs(
        gradients[0] + 1j*gradients[1]).real, axis=(1, 2)) # [TD]
    return tenengrad.mean()

class SSIM():
    def __init__(self, window_size=7, k1=0.01, k2=0.03, dim=2):
        assert window_size % 2
        self.window_size = window_size
        self.k1 = k1
        self.k2 = k2
        self.dim = dim
        assert dim == 2 or dim == 4

    def __call__(self, target, recon, alpha=1, beta=1, gamma=1): # [T, D, H, W, 2] or [TD, H, W, 2]
        target = torch.abs(torch.view_as_complex(target))
        recon = torch.abs(torch.view_as_complex(recon))
        
        if len(target.shape) == 4:
            nbin = target.shape[0]
            slice_dims = target.shape[2:]
        elif len(target.shape) == 3:
            nbin = 1
            slice_dims = target.shape[1:]

        exp_t, exp_r, var_t, var_r, covar, c1, c2 = self.compute_statistics(
            target.reshape(-1, 1, *slice_dims),
            recon.reshape(-1, 1, *slice_dims), nbin=nbin
        )

        if alpha == 1 and beta == 1 and gamma == 1:
            ssim = self.compute_ssim(
                exp_t, exp_r, var_t, var_r, covar, c1, c2)
            return ssim[~ssim.isnan()].mean()
        
        luminosity = self.compute_luminosity(exp_t, exp_r, c1)
        contrast, structure = self.compute_contrast_structure(
            var_t, var_r, covar, c2)
        ssim = luminosity ** alpha * contrast ** beta * structure ** gamma
        return ssim.mean()

    def compute_luminosity(self, exp_t, exp_r, c1):
        return (
            2 * exp_t * exp_r + c1) / (exp_t ** 2 + exp_r ** 2 + c1)

    def compute_contrast_structure(self, var_t, var_r, covar, c2): # [TD, 1, h, w] for 2d case, c2: [TD, 1, 1, 1]
        contrast = (2 * var_t * var_r + c2) / (var_t + var_r + c2)
        structure = (covar + c2 / 2 ) / (torch.sqrt(var_t * var_r) + c2 / 2)
        return contrast, structure # [TD, 1, h, w]

    def compute_ssim(self, exp_t, exp_r, var_t, var_r, covar, c1, c2):
        return (2 * exp_t * exp_r + c1) * (2 * covar + c2) / (
            (exp_t ** 2 + exp_r ** 2 + c1) * (var_t + var_r + c2)
        )
    
    def compute_statistics(
        self, target, recon, window_size=None, nbin=None): # [TD, 1, H, W]
        window_size = self.window_size if window_size is None else window_size
        
        cov_norm = window_size ** self.dim / (window_size ** self.dim - 1)
        
        exp_t = torch.nn.functional.avg_pool2d(target, window_size, stride=1) # [TD, 1, h, w]
        exp_r = torch.nn.functional.avg_pool2d(recon, window_size, stride=1)
        exp_tr = torch.nn.functional.avg_pool2d(
            target * recon, window_size, stride=1)
        exp_tt = torch.nn.functional.avg_pool2d(
            target * target, window_size, stride=1)
        exp_rr = torch.nn.functional.avg_pool2d(
            recon * recon, window_size, stride=1)

        if self.dim == 4:
            exp_t = exp_t.reshape(
                nbin, int(target.shape[0] / nbin), *exp_t.shape[2:] # [T, D, h, w] or [1, TD, h, w]
            ).moveaxis((2, 3), (0, 1)) # [h, w, T, D] or [h, w, 1, TD]
            exp_r = exp_r.reshape(
                nbin, int(recon.shape[0] / nbin), *exp_r.shape[2:] 
            ).moveaxis((2, 3), (0, 1))
            exp_tr = exp_tr.reshape(
                nbin, int(recon.shape[0] / nbin), *exp_tr.shape[2:] 
            ).moveaxis((2, 3), (0, 1))
            exp_tt = exp_tt.reshape(
                nbin, int(recon.shape[0] / nbin), *exp_tt.shape[2:] 
            ).moveaxis((2, 3), (0, 1))
            exp_rr = exp_rr.reshape(
                nbin, int(recon.shape[0] / nbin), *exp_rr.shape[2:] 
            ).moveaxis((2, 3), (0, 1))
            if nbin == 1:
                avgpool = torch.nn.functional.avg_pool1d
                shape = (-1, *exp_t.shape[2:])
            else:
                avgpool = torch.nn.functional.avg_pool2d
                shape = (-1, 1, *exp_t.shape[2:])
            exp_t = exp_t.reshape(*shape)
            exp_r = exp_r.reshape(*shape)
            exp_tr = exp_tr.reshape(*shape)
            exp_tt = exp_tt.reshape(*shape)
            exp_rr = exp_rr.reshape(*shape) # [hw, 1, T, D] or [hw, 1, TD]
            
            exp_t = avgpool(exp_t, window_size, stride=1) # [hw, 1, t, d] or [hw, 1, td]
            exp_r = avgpool(exp_r, window_size, stride=1)
            exp_tr = avgpool(exp_tr, window_size, stride=1)
            exp_tt = avgpool(exp_tt, window_size, stride=1)
            exp_rr = avgpool(exp_rr, window_size, stride=1)

        var_t = cov_norm * (exp_tt - exp_t * exp_t)
        var_r = cov_norm * (exp_rr - exp_r * exp_r)
        covar = cov_norm * (exp_tr - exp_t * exp_r)

        if self.dim == 2:
            dynamic_range = target.reshape(target.shape[0], -1).max(axis=1)[0] # [TD]
            dynamic_range = dynamic_range.reshape(
                *dynamic_range.shape, 1, 1, 1) # [TD, 1, 1, 1]
        else:
            dynamic_range = target.max()
        c1 = (self.k1 * dynamic_range) ** 2
        c2 = (self.k2 * dynamic_range) ** 2

        return exp_t, exp_r, var_t, var_r, covar, c1, c2

class MS_SSIM(SSIM):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        assert self.dim == 2

    def __call__(self, target, recon): # [T, D, H, W, 2] or [D, H, W, 2]
        is_with_dynamic_gaps = len(target.shape) == 5
        exponents = [0.0448, 0.2856, 0.3001, 0.2363, 0.1333]
        target = torch.abs(torch.view_as_complex(target))
        recon = torch.abs(torch.view_as_complex(recon))

        if len(target.shape) == 4:
            slice_dims = target.shape[2:]
        else:
            slice_dims = target.shape[1:]
        t = target.reshape(-1, 1, *slice_dims) # [TD, 1, H, W]
        r = recon.reshape(-1, 1, *slice_dims)
        contrast_product, structure_product = 1, 1
        
        for exponent in exponents[:-1]:
            _, _, var_t, var_r, covar, _, c2 = self.compute_statistics(t, r)

            contrast, structure = self.compute_contrast_structure(
                var_t, var_r, covar, c2
            ) # [TD, 1, h, w]
            
            contrast_product = contrast_product * torch.pow(
                contrast.mean(dim=(1, 2, 3)), exponent)
            structure_product = structure_product * torch.pow(
                structure.mean(dim=(1, 2, 3)), exponent)

            t = torch.nn.functional.avg_pool2d(t, 2) # [TD, 1, h, w]
            r = torch.nn.functional.avg_pool2d(r, 2)
        
        if t.shape[-1] < self.window_size or t.shape[-2] < self.window_size:
            window_size = min(t.shape[-1], t.shape[-2])
        else:
            window_size = self.window_size
        exp_t, exp_r, var_t, var_r, covar, c1, c2 = self.compute_statistics(
            t, r, window_size=window_size)
        contrast, structure = self.compute_contrast_structure(
            var_t, var_r, covar, c2)
        contrast_product = contrast_product * torch.pow(
            contrast.mean(dim=(1, 2, 3)), exponents[-1])
        structure_product = structure_product * torch.pow(
            structure.mean(dim=(1, 2, 3)), exponents[-1])
        luminosity = torch.pow(
            self.compute_luminosity(exp_t, exp_r, c1).mean(dim=(1, 2, 3)),
            exponents[-1]
        )
        
        ms_ssim = luminosity * contrast_product * structure_product
        return ms_ssim.mean()

def ms_ssim(y, x): # [T, D, H, W, 2] or [D, H, W, 2]
    x = torch.abs(torch.view_as_complex(x))
    y = torch.abs(torch.view_as_complex(y))
    kernel_size = 11
    if x.shape[-1] <= 160 or x.shape[-2] <= 160:
        kernel_size = 9
    metric = torchmetrics.image.MultiScaleStructuralSimilarityIndexMeasure(
        gaussian_kernel=False, kernel_size=kernel_size).to(device=y.device)
    
    return metric(
        x.reshape(-1, 1, *x.shape[-2:]),
        y.reshape(-1, 1, *y.shape[-2:])
    )

def plot_images(
    config, target, recon, dataset, maskname, checkpoint, writer, other_gif_name = None):
    
    def transform_image_pair(t, x):
        minval = min(t.min().item(), x.min().item(), 0)
        shifted_t = t - minval
        shifted_x = x - minval
        shifted_max = max(shifted_t.max().item(), shifted_x.max().item())
        return shifted_t / shifted_max, shifted_x / shifted_max
    def make_false_color_components(t, x):
        transformed_t, transformed_x = transform_image_pair(t, x)
        return torch.stack((transformed_t, transformed_x, transformed_t), 0)
    def normalize(x):
        return x / x.max()

    nrow = 5 if not config['volume-slices'] else len(
        config['volume-slices'].split())
    if not config['volume-slices']:
        vol_slices = np.linspace(0, recon.size(1) - 1, num=nrow).astype(int)
    else:
        vol_slices = [int(vs) for vs in config['volume-slices'].split()]
    if not config['bin-slices']:
        bin_slices = np.linspace(0, 9, num=nrow).astype(int)
    else:
        bin_slices = [int(bs) for bs in config['bin-slices'].split()]
        
    #targets_width_by_height = [torch.rot90(
    #    target[i, j].abs(), 3, [0, 1]) for i, j in zip(bin_slices, vol_slices)]
    recons_width_by_height = [torch.rot90(
        recon[i, j].abs(), 3, [0, 1]) for i, j in zip(bin_slices, vol_slices)]

    #writer.add_image(f'Width-by-Height-errors/{dataset}/{maskname}', torchvision.utils.make_grid(torch.stack([make_false_color_components(t, r) for t, r in zip(targets_width_by_height, recons_width_by_height)], 0), nrow=nrow), int(checkpoint))

    width_by_height = torchvision.utils.make_grid(torch.stack([normalize(r) for r in recons_width_by_height], 0).unsqueeze(1), nrow=nrow) #[normalize(t) for t in targets_width_by_height] +
    writer.add_image(f'Width-by-Height/{dataset}/{maskname}', width_by_height, int(checkpoint))
    #recon = torch.stack([torch.rot90(recon[:, i, :, :], 3, [1, 2]) for i in range(recon.size(1))], dim=1)
    recon = torch.abs(recon)
    recon = (recon - recon.min()) / (recon.max() - recon.min())
    recon = recon.unsqueeze(0)
    recon = recon.repeat(1, 1, 3, 1, 1)

    gif_name = 'gif'
    if other_gif_name != None:
        gif_name = other_gif_name
    writer.add_video(gif_name, recon, int(checkpoint), fps = 30)
    writer.flush()
    #writer.close()