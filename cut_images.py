
'''
Test files:

/home/rw264/Documents/SDSS/Sunayana_Files/fits/frame-g-001462-2-0371.fits.bz2
/home/rw264/Documents/SDSS/Sunayana_Files/fits/frame-i-001462-2-0371.fits.bz2
/home/rw264/Documents/SDSS/Sunayana_Files/fits/frame-r-001462-2-0371.fits.bz2
/home/rw264/Documents/SDSS/Sunayana_Files/fits/frame-u-001462-2-0371.fits.bz2
/home/rw264/Documents/SDSS/Sunayana_Files/fits/frame-z-001462-2-0371.fits.bz2
'''

import sys


import matplotlib
matplotlib.use('agg')
import aplpy
import pyregion as pr
from astropy import units as u
from astropy import coordinates,wcs
from astropy.nddata import Cutout2D
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import math as m
