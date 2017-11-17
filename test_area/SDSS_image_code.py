import os
import pyregion as pr
from astropy import units as u
from astropy import coordinates,wcs
from astropy.nddata import Cutout2D
from astropy.io import fits
import aplpy
import matplotlib.pyplot as plt
import numpy as np
import math as m

def get_ra_dec(coords):
    with open(coords) as f:
        positions=f.readlines()
    return positions

def get_votable(ra,dec,acc,bands,work,outfile):
    os.chdir(work)
    link = "http://skyserver.sdss.org/dr14/SkyServerWS/SIAP/getAllSIAPInfo?POS={},{}&SIZE={}&FORMAT=image/fits&bandpass={}".format(ra,dec,acc,bands)
    os.system("wget '{}' -O {}".format(link,outfile))
    #print "wget '{}' -O {}".format(link,"temp.xml")

def get_link_from_vot(infile,outfile):
    os.system("grep -r '<url>http' {} > temp.txt".format(infile))
    os.system("grep -r 'frame-g\|frame-r\|frame-i' temp.txt > {}".format(infile))
    os.system("sed 's/<\/url>//' {} > temp.txt".format(infile))
    os.system("sed 's/<url>//' temp.txt >{} ".format(outfile))
    os.system("rm -rf temp.txt")



def read_links(links):
    with open(links) as f:
        contents=f.readlines()

    for line in range(0,len(contents)):
        temp = contents[line]
        contents[line] = temp.replace("\r\n","").replace(" ","")
    return contents



def get_files(links,ra,dec,output_dir):
    output_filenames={}
    try:
        os.chdir(output_dir)
    except OSError:
        os.mkdir(output_dir)
        os.chdir(output_dir)

    for alink in links:
        filename=
        os.system("wget {}".format(alink))




def pull_coords(output_filenames):
    ralong,declong = output_filenames.keys()[0].split("_")
    ra = ralong.split("=")[1]
    dec = declong.split("=")[1]
    return ra,dec

def make_cutout(ra,dec,size_arc,filename):
    fit_image = fits.open(filename)
    header = fit_image[0].header
    w = wcs.WCS(fit_image[0].header)
    im1 = fit_image[0].data
    center = coordinates.SkyCoord(ra,dec,unit='deg')
    size = u.Quantity([size_arc,size_arc],u.arcmin)
    co = Cutout2D(im1,center,size,wcs=w)
    return co,header

def plot_fits(data,header,filename):
    hdu = fits.PrimaryHDU(data=data,header=header)
    gc=aplpy.FITSFigure(hdu)
    gc.show_grayscale(vmin=0, vmax=5.0)
    gc.hide_axis_labels()
    gc.ticks.hide()
    gc.hide_tick_labels()
    gc.frame.set_linewidth(0)
    gc.save(filename)

def save_fits(data,header,filename):
    hdu = fits.PrimaryHDU(data=data,header=header)
    fits.writeto(filename,data=data,header=header)

if __name__ == '__main__':
    links = "/home/rw264/Documents/SDSS/test_area/link.txt"
    cords = "/home/rw264/Documents/SDSS/test_area/coords.txt"
    storage = '/home/rw264/Documents/SDSS/test_area/'
    xml = "temp.xml"
    linkfile = "link.txt"
    positions = get_ra_dec(cords)
    ra,dec = positions[0].replace("\n","").split(",")
    get_votable(ra,dec,0.1,"gri",storage,xml)
    get_link_from_vot(xml,linkfile)
    links = read_links(linkfile)
    links
    ra,dec = pull_coords(filenames)
    fits_file = filenames[filenames.keys()[0]][4]
    try:
        cut_out,c0_header = make_cutout(ra,dec,12,fits_file)
    #plot_fits(cut_out.data,c0_header,"test_of_code_z.png")
        save_fits(cut_out.data,c0_header,"cutout_"+fits_file[0:-4])
    except ValueError:
        None
