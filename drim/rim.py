__author__ = 'Kai'

import torch
from torch import nn
import numpy as np

import logging

def set_conv_block(module, name, n_chan, n_feat, dilation):
    setattr(module, name, nn.Sequential())
    n_layer = -1
    for kernel in module.kernel:
        if not kernel is None:
            n_layer += 1
    for i, kernel in enumerate(module.kernel):
        if not kernel is None:
            if kernel[0] > 1:
                padt = int((
                    kernel[0] + (kernel[0] - 1)*(dilation - 1)) / 2)
                getattr(module, name).append(CyclicPad1d(padt))
            if kernel[1] > 1 or kernel[2] > 1:
                padw = int((
                    kernel[2] + (kernel[2] - 1)*(dilation - 1)) / 2)
                padd = int((
                    kernel[1] + (kernel[1] - 1)*(dilation - 1)) / 2)
                getattr(module, name).append(
                    nn.ReplicationPad3d((padw, padw, padd, padd, 0, 0))
                )
            getattr(module, name).append(nn.Conv3d(
                n_chan, n_feat, kernel, dilation=dilation, bias=True))
            if not (name == 'conv_out' and i == n_layer):
                getattr(module, name).append(nn.ReLU(inplace=True))
            n_chan = n_feat

class RecurrentInferenceMachine(nn.Module):
    def __init__(
        self, nfeature=128, kernel=((1, 3, 3), (3, 1, 1), None),
        temporal_rnn=False, mute=True):
        super().__init__()
        self.logger = logging.getLogger(type(self).__name__)
        self.nfeature = nfeature
        if kernel[2] is not None: assert kernel[1] is not None
        if kernel[1] is not None: assert kernel[0] is not None
        self.kernel = kernel
        set_conv_block(self, 'conv_in', 4, self.nfeature, 1)
        self.recurrent1 = LatticeRecurrentUnit(
            self.nfeature, self.nfeature, temporal=temporal_rnn, mute=mute)
        set_conv_block(self, 'conv_between', self.nfeature, self.nfeature, 1)
        self.recurrent2 = LatticeRecurrentUnit(
            self.nfeature, self.nfeature, temporal=temporal_rnn, mute=mute)
        set_conv_block(self, 'conv_out', self.nfeature, 2, 1)
        if not mute:
            self.logger.info("Initialized a recurrent inference machine of "
                  f"{get_num_params(self)} parameters. Convolutional blocks "
                  f"are conv_in: {self.conv_in}, conv_between: "
                  f"{self.conv_between}, and conv_out: {self.conv_out}.")

    def forward(self, x, hidden): # [B, F, T, D, W], [B, F, T, D, W]
        x = self.conv_in(x)
        x, hidden[0] = self.recurrent1(x, hidden[0])
        x = self.conv_between(x)
        x, hidden[1] = self.recurrent2(x, hidden[1])
        x = self.conv_out(x)
        return x, hidden # [B, T, D, W], [B, F, T, D, W]

class GradRim(nn.Module):
    def __init__(self, fourier_dim=-1):
        super().__init__()
        self.fourier_dim = fourier_dim

    def forward(self, x, measurements, sense, mask): # [B, T, D, W, 2], 3 * [B, C, T, D, W, 2]
        x = x.unsqueeze(1)
        x = torch.view_as_complex(torch.stack((
            sense[..., 0] * x[..., 0] - sense[..., 1] * x[..., 1],
            sense[..., 0] * x[..., 1] + sense[..., 1] * x[..., 0]), -1))
        kspace = torch.view_as_real(torch.fft.fftn(x, dim=self.fourier_dim))
        error = torch.where(
            mask,
            kspace - measurements,
            torch.Tensor([0]).to(measurements))
        error = torch.view_as_real(
            torch.fft.ifftn(
                torch.view_as_complex(error), dim=self.fourier_dim))
        gradient = torch.stack((
            torch.sum(
                sense[..., 0] * error[..., 0] + sense[..., 1] * error[..., 1],
                1),
            torch.sum(
                sense[..., 0] * error[..., 1] - sense[..., 1] * error[..., 0],
                1)
            ), 1)
        return gradient # [B, T, D, W]

class CyclicPad1d(nn.Module):
    def __init__(self, n):
        super().__init__()
        self.n = n
    def __repr__(self):
        return f"{type(self).__name__}: ({self.n}, {self.n}, 0, 0, 0, 0)"
    def forward(self, input): # [B, F, T, D, W]
        return torch.cat(
            (input[:, :, -self.n:], input, input[:, :, :self.n]), 2)

class LatticeRecurrentUnit(nn.Module):
    def __init__(self, input_size, hidden_size, temporal=False, mute=True):
        super().__init__()
        self.logger = logging.getLogger(type(self).__name__)
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.temporal = temporal
        self.conv_ih = nn.Conv3d(
            input_size, 6 * hidden_size,
            (3 if temporal else 1, 1, 1), bias=True)
        self.conv_hh = nn.Conv3d(
            hidden_size, 6 * hidden_size,
            (3 if temporal else 1, 1, 1), bias=False)
        self.bias_h = nn.parameter.Parameter(
            torch.zeros(self.hidden_size).view(1, self.hidden_size, 1, 1, 1))
        if not mute:
            self.logger.info(
                f'{type(self).__name__} '
                f'{"with temporal convolutions" if temporal else ""}'
                f'initialized with input size {input_size}, '
                f'hidden size {hidden_size}, containing '
                f'{get_num_params(self)} parameters.')

    def forward(self, input, hidden):
        if self.temporal:
            padded_input = CyclicPad1d(1)(input)
            padded_hidden = CyclicPad1d(1)(hidden)
        else:
            padded_input = input
            padded_hidden = hidden
        ii_r, ih_r, ii_z, ih_z, ii_c, ih_c = self.conv_ih(
            padded_input).chunk(6, 1)
        hi_r, hh_r, hi_z, hh_z, hi_c, hh_c = self.conv_hh(
            padded_hidden).chunk(6, 1)
        reset_i = torch.sigmoid(ii_r + hi_r)
        reset_h = torch.sigmoid(ih_r + hh_r)
        update_i = torch.sigmoid(ii_z + hi_z)
        update_h = torch.sigmoid(ih_z + hh_z)
        candidate_i = torch.tanh(ih_c + reset_h * (hh_c + self.bias_h))
        candidate_h = torch.tanh(reset_i * ii_c + hi_c)
        return (update_i * candidate_h + (1 - update_i) * input,
                update_h * candidate_i + (1 - update_h) * hidden)

class InitRim(nn.Module):
    def __init__(
        self, x_ch, out_channels, kernel, multiscale_depth=0, mute=True):
        """
        Learned initializer for RIM, based on multi-scale context aggregation 
        with dilated convolutions, that replaces zero initializer for the RIM
        hidden vector. Inspired by "Multi-Scale Context Aggregation by
        Dilated Convolutions" (https://arxiv.org/abs/1511.07122)
        
        Parameters
        ----------
        x_ch : int
            Input channels.
        out_channels : List of int
            Number of features for the hidden states in the RIM.
        dilations: tuple
            Dilations of the convolutional layers of the initializer.
        multiscale_depth : 1
            Number of feature layers to aggregate for the output, if 1,
            multi-scale context aggregation is disabled.
        """
        super().__init__()
        self.logger = logging.getLogger(type(self).__name__)
        self.conv_blocks = nn.ModuleList()
        self.out_blocks = nn.ModuleList()
        self.multiscale_depth = multiscale_depth
        self.kernel = kernel
        tch = x_ch
        outchannels = [
            max(out_channels)//4] + 2*[
            max(out_channels)//2] + [
            max(out_channels)]
        for nconv, (c, d) in enumerate(
            zip(outchannels, [1] + [2**r for r in range(3)])):
            set_conv_block(self, f'conv{nconv+1}', tch, c, d)
            self.conv_blocks.append(getattr(self, f'conv{nconv+1}'))
            tch = c
        tch = np.sum(outchannels[-multiscale_depth:])
        for ch in out_channels:
            outblock = nn.Sequential(
                nn.Conv3d(tch, ch, 1),
                nn.ReLU(inplace=True))
            self.out_blocks.append(outblock)
        if not mute: 
            self.logger.info(
                f'{type(self).__name__} initialized, containing '
                f'{get_num_params(self)} parameters.')

    def forward(self, x): # [B, F, T, D, W]
        features = []
        for block in self.conv_blocks:
            x = block(x)
            if self.multiscale_depth != 1:
                features.append(x)
        if self.multiscale_depth != 1:
            x = torch.cat(features[-self.multiscale_depth:], dim=1)
        outlst = []
        for block in self.out_blocks:
            outlst.append(block(x))
        return list(outlst) # 2x [B, F, T, D, W]

def get_num_params(module):
    n = 0
    for na, p in module.named_parameters():
        module.logger.debug(f"{na}: {p.size()}")
        n += p.numel()
    return n