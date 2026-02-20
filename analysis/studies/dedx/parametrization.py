import os
import sys
import json
import numpy as np


# load parameters
with open('params.json', 'r') as f:
    PARAMS = json.load(f)


def get_parametrized_curve(p, species='proton', subsystem='wires'):
    params = PARAMS[species][subsystem]
    a, b, c, d, e = params
    curve = a * (b + c*np.log(p) + d*np.power(p, e))
    return curve
