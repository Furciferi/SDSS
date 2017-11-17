import sys
########### Reese's library ############
sys.path.append("/home/r/rw/rw264/Workspace/lib/python2.7/site-packages")

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


fit_image = fits.open('0674380301.fits')
w = wcs.WCS(fit_image[0].header)
im1 = fit_image[0].data


region_file = "0674380301.reg"

r = pr.open(region_file)
for i in range(0,1):#len(r)):
xapa_psf = r[i].coord_list
if xapa_psf[2]!=xapa_psf[3]:
continue


center = coordinates.SkyCoord.from_pixel(xapa_psf[0],xapa_psf[1],w,origin=0,mode='wcs')
size = coordinates.SkyCoord.from_pixel(xapa_psf[0]+xapa_psf[2],xapa_psf[1],w,origin=0,mode='wcs')
size = center.separation(size).arcsec
co = Cutout2D(im1,center,4*size*u.arcsec,wcs=w)

hdu = fits.PrimaryHDU(data=co.data,header=co.wcs.to_header())
gc=aplpy.FITSFigure(hdu)
gc.show_grayscale(vmin=0, vmax=5.0)
gc.hide_axis_labels()
gc.hide_tick_labels()
gc.ticks.hide()
gc.frame.set_linewidth(0)
gc.show_regions(region_file)
gc.save('0674380301_psf_cutout'+str(i)+'.png')
