import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame as pdf
import sys


def make_quality_plot(ifname=''):

    # read .csv file
    df = pd.read_csv(ifname, comment='#')
