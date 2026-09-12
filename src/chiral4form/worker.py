"""Run one configured prime as a scheduler-safe isolated pipeline shard.

Each prime receives a separate output/checkpoint directory. No implicit
background process is started and no upstream source code is executed.
"""
from copy import deepcopy
from pathlib import Path
from .pipeline import run,validate_config


def run_prime(config,prime_index,root,output,source=None,progress=print):
    validate_config(config)
    if type(prime_index) is not int or not 0<=prime_index<len(config['primes']):
        raise ValueError('prime-index outside the configured list')
    one=deepcopy(config);one['primes']=[config['primes'][prime_index]]
    return run(one,root,Path(output)/f'prime-{one["primes"][0]}',source,progress)
