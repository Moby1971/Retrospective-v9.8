# encoding: utf-8
__author__ = 'Jonas'

import logging
import sys

def setup(
    use_stdout=True, filename=None, log_level=logging.INFO,
    redirect_stdout=False, redirect_stderr=False):
    """setup some basic logging"""

    log = logging.getLogger('')
    log.setLevel(log_level)
    fmt = logging.Formatter(
        "%(levelname)-5s - %(asctime)s [%(name)-15s] %(message)s",
        datefmt="%y-%m-%d %H:%M:%S")

    if use_stdout:
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(log_level)
        ch.setFormatter(fmt)
        log.addHandler(ch)

    if not filename is None:
        fh = logging.FileHandler(filename)
        fh.setLevel(log_level)
        fh.setFormatter(fmt)
        log.addHandler(fh)

    if redirect_stderr:
        stderr_logger = logging.getLogger('STDERR')
        sl = StreamToLogger(stderr_logger, logging.ERROR)
        sys.stderr = sl

    if redirect_stdout:
        stdout_logger = logging.getLogger('STDOUT')
        sl = StreamToLogger(stdout_logger, logging.ERROR)
        sys.stdout = sl  

class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass