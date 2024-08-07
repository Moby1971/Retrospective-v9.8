__author__ = 'Kai'

import sys
sys.path.insert(0, './drim')
import os
import configparser
import logging
from logger import setup as log_setup
from train import train_model
from validate import validate_model
from timer import time_model
from reconstruct import reconstruct

# sys.argv[1] = reconstruct
# sys.argv[2] = train directory
# sys.argv[3] = temporary directory where mat files are located
# sys.argv[4] = checkpoint 

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if '-h' in sys.argv or '--help' in sys.argv:
        print(
            "Script contains four different programs:"
            "\t1) train: train a network."
            "\t2) validate: validate a network."
            "\t3) time: time a network."
            "\t4) reconstruct: reconstruct data using a pre-trained network."
            "For more details, specify which program you need help with.")
    config = configparser.ConfigParser()
    config.read('drim/configfile-example')
    if sys.argv[1] == "train" or sys.argv[1] == 'validate' or (
        sys.argv[1] == 'reconstruct' or sys.argv[1] == 'time'):
        for arg in sys.argv:
            if arg.startswith('@'):
                config.read(arg[1:])
                break

    if sys.argv[1] == 'train':
        if 'saved-model' in config['train']:
            traindir = os.path.dirname(os.path.dirname(
                config['train']['saved-model']))
            new_config_name = config['train']['saved-model']
            new_config_name = new_config_name.split('checkpoint')[1][:-3] 
            new_config_name = new_config_name + '.ini'
            with open(os.path.join(traindir, new_config_name), 'w') as f:
                config.write(f)
        else:
            traindir = config['train']['train-dir']
            if os.path.exists(traindir):
                num = 2
                traindir = config['train']['train-dir'] + f'-{num}'
                while os.path.exists(traindir):
                    num += 1
                    traindir = config['train']['train-dir'] + f'-{num}'
            config['train']['train-dir'] = traindir
            os.makedirs(os.path.join(traindir, 'logs'))
            os.mkdir(os.path.join(traindir, 'network-parameters'))
            with open(
                os.path.join(config['train']['train-dir'], 'config.ini'), 'w'
            ) as f:
                config.write(f)
        
        log_setup(
            True,
            os.path.join(traindir, 'logs', 'train_run_log.txt'),
            log_level=logging.INFO)
        logger.info(f"Storing train run in {config['train']['train-dir']}")
        for key, val in config['train'].items():
            logger.info(f"{key}: {val}")

        train_model(config['train'])
        logger.info("Finished training!")

    if sys.argv[1] == 'validate':
        log_setup(
            True,
            os.path.join(config['validate']['train-dir'],
                         'logs', 'validation_log.txt'),
            log_level=logging.INFO)
        for key, val in config['validate'].items():
            logger.info(f"{key}: {val}")
        validate_model(config['validate'])
        logger.info('Finished validation.')

    if sys.argv[1] == 'time':
        if os.path.exists(config['time']['train-dir']):
            time_log = os.path.join(
                config['time']['train-dir'], 'logs', 'time_log.txt')
        else:
            if not os.path.exists('time_logs'):
                os.mkdir('time_logs')
            time_log = os.path.join('time_logs', config['time']['train-dir'])
        log_setup(
            True,
            time_log,
            log_level=logging.INFO)
        for key, val in config['time'].items():
            logger.info(f"{key}: {val}")
        time_model(config['time'])
        logger.info('Finished timing.')

    if sys.argv[1] == 'reconstruct':
        if '-h' in sys.argv or '--help' in sys.argv:
            pass
        config = config['reconstruct']
        if os.path.exists(config['wdir']):
            raise ValueError(f'wdir {config["wdir"]} already exists.')
        else:
            log_setup(
                True,
                os.path.join(sys.argv[2],
                         'logs', 'reconstruction_log.txt'),
                log_level=logging.INFO)
            logger.info('Reconstructing...')
            reconstruct(config)