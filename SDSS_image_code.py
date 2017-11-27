import os
import pyregion as pr
from astropy import units as u
from astropy import coordinates,wcs
from astropy.nddata import Cutout2D,utils
from astropy.io import fits
import aplpy
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.ndimage

def get_ra_dec(coords):
    with open(coords) as f:
        positions=f.readlines()
    return positions

def get_votable(ra,dec,acc,bands,work,outfile):
    os.chdir(work)
    link = "http://skyserver.sdss.org/dr14/SkyServerWS/SIAP/getSIAPInfo?POS={},{}&SIZE={}&FORMAT=image/fits&bandpass={}".format(ra,dec,acc,bands)
    print "\n {} \n".format(link)
    os.system("aria2c -x 4 --allow-overwrite '{}' -o {}".format(link,outfile))


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



def get_files(links,ra,dec,output_dir,work,download):
    dir_name="/ra_{}_dec_{}".format(ra,dec)
    output_filenames=[]
    try:
        os.chdir(output_dir+dir_name)
    except OSError:
        os.mkdir(output_dir+dir_name)
        os.chdir(output_dir+dir_name)

    for alink in links:
        if download==True:
            os.system("aria2c -x 4 {}".format(alink))
        output_filenames.append(alink.split("/")[-1])
    os.chdir(work)
    return output_filenames



def pull_coords(output_filenames):
    ralong,declong = output_filenames.keys()[0].split("_")
    ra = ralong.split("=")[1]
    dec = declong.split("=")[1]
    return ra,dec

def make_cutout(ra,dec,size_arc,filename,storage,_mode):
    os.chdir(storage+"/ra_{}_dec_{}".format(ra,dec))
    hdulist = fits.open(filename)
    header = hdulist[0].header
    header["RADESYS"]=header["RADECSYS"]
    header.pop("RADECSYS")
    #print header["CRVAL1"],header["CRVAL2"]
    w = wcs.WCS(hdulist[0].header)
    im1 = hdulist[0].data
    center = coordinates.SkyCoord(ra,dec,unit='deg')
    #print center
    size = u.Quantity([size_arc,size_arc],u.pixel)
    if _mode == "partial":
        co =Cutout2D(im1,center,size,wcs=w,mode='partial',fill_value='Nan')
    else:
        co =Cutout2D(im1,center,size,wcs=w,mode='strict',fill_value='Nan')
    print co.center_cutout
    header["CRVAL1"]=ra
    header["CRVAL2"]=dec
    header["CRPIX1"]=size_arc//2
    header["CRPIX2"]=size_arc//2
    return co,header

def plot_fits(data,header,filename,work,storage):
    os.chdir(storage+"/ra_{}_dec_{}".format(ra,dec))
    hdu = fits.PrimaryHDU(data=data,header=header)
    gc=aplpy.FITSFigure(hdu)
    gc.show_grayscale(vmin=0, vmax=5.0)
    gc.hide_axis_labels()
    gc.ticks.hide()
    gc.hide_tick_labels()
    gc.frame.set_linewidth(0)
    os.chdir(work)
    gc.save(filename)

def save_fits(data,header,filename,output):
    try:
        os.chdir(output)
    except OSError:
        os.mkdir(output)
        os.chdir(output)
    hdu = fits.PrimaryHDU(data=data,header=header)
    fits.writeto(filename,data=data,header=header)

def make_and_save_cuts(filenames,ra,dec,output,size,storage,skipto,_mode):
    if skipto!=[]:
        filenames=skipto
        skipto=[]
    for name in filenames:
        print str(100*round(float(filenames.index(name))/float(len(filenames)),4))+"%"
        try:
            cut_out,header = make_cutout(ra,dec,size,name,storage,_mode)
            save_fits(cut_out.data,header,name[:-4],output)
            skipto.append(name)
        except (utils.NoOverlapError,IOError,utils.PartialOverlapError):

            continue
    return skipto

def swarp_fits(output,size):
    for band in ["*frame-g*","*frame-i*","*frame-r*"]:
        os.chdir(output)
        os.system("swarp {} -IMAGE_SIZE '{},{}' -IMAGEOUT_NAME '{}.fits'".format(band,size,size,band[-2:-1]))
        os.system("rm swarp.xml")


def make_the_coads(filenames,ra,dec,output,storage):
    size=455
    try:
        _3x3=make_and_save_cuts(filenames,ra,dec,output+"_{}-arcmin/".format(size),size,storage,filenames,'strict')
        swarp_fits(output+"_{}-arcmin".format(size),size)
    except (OSError,IOError):
        print "There was an error loading \nimage_size={} \n".format(size)
    size=909
    try:
        _6x6=make_and_save_cuts(filenames,ra,dec,output+"_{}-arcmin/".format(size),size,storage,_3x3,'strict')
        swarp_fits(output+"_{}-arcmin".format(size),size)
    except (OSError,IOError):
        print "There was an error loading \nimage_size={} \n".format(size)
    size=1365
    try:
        _9x9=make_and_save_cuts(filenames,ra,dec,output+"_{}-arcmin".format(size),size,storage,_6x6,'strict')
        swarp_fits(output+"_{}-arcmin".format(size),size)
    except (IOError,OSError):
        print "There was an error loading \nimage_size={} \n Now creating an Image with filled in Blank space...".format(size)
        try:
            make_and_save_cuts(filenames,ra,dec,output+"_{}-arcmin".format(size),size,storage,_6x6,'partial')
            swarp_fits(output+"_{}-arcmin".format(size),size)
        except (IOError,OSError):
            print "Partial cutout failed... \n\n"
    size=1818
    try:
        _12x12=make_and_save_cuts(filenames,ra,dec,output+"_{}-arcmin".format(size),size,storage,_6x6,'strict')
        swarp_fits(output+"_{}-arcmin".format(size),size)
    except (IOError,OSError):
        print "There was an error loading \nimage_size={} \n Now creating an Image with filled in Blank space...".format(size)
        try:
            make_and_save_cuts(filenames,ra,dec,output+"_{}-arcmin".format(size),size,storage,_6x6,'partial')
            swarp_fits(output+"_{}-arcmin".format(size),size)
        except (IOError,OSError):
            print "Partial cutout failed... \n\n"

def final_images(xcs_coords,peak_coords,output,size,fileno):
    try:
        RA = float(xcs_coords[fileno].replace("\n","").split("\t")[0])
        DEC = float(xcs_coords[fileno].replace("\n","").split("\t")[1])
        RA_2= float(peak_coords[fileno].replace("\n","").split("\t")[0])
        DEC_2= float(peak_coords[fileno].replace("\n","").split("\t")[1])
        RA_3=float(ra)
        DEC_3=float(dec)
        os.chdir(output+"_{}-arcmin".format(size))
        os.system("stiff i.fits r.fits g.fits")
        oi =aplpy.FITSFigure("i.fits",dimensions=[0,1],slices=[0])
        oi.show_rgb("stiff.tif")
        oi.axis_labels.hide()
        oi.tick_labels.hide()
        oi.ticks.hide()
        oi.show_markers(RA, DEC, edgecolor='cyan',linestyle='--',marker="o",s=5000)
        oi.show_markers(RA_2, DEC_2, edgecolor='magenta',linestyle='--',marker="o",s=5000)
        oi.show_markers(RA_3, DEC_3, edgecolor='yellow',linestyle='--',marker="o",s=5000)
        oi.save("final-{}-arcmin.jpeg".format(size),dpi=100)
    except (OSError,IOError):
        print "There was an error loading \nimage_size={} \n".format(size)

if __name__ == '__main__':
    links = "/home/rw264/Documents/SDSS/test_area/link.txt"
    cords = "/home/rw264/Documents/SDSS/test_area/coords.txt"
    work = '/home/rw264/Documents/SDSS/test_area/'
    storage = '/home/rw264/Documents/Sunayana_Files/fits'

    xml = "temp.xml"
    linkfile = "link.txt"
    positions = get_ra_dec(cords)
    fileno = 2
    ra,dec = positions[fileno].replace("\n","").split(",")
    get_votable(ra,dec,"0.1","gri",work,xml)
    get_link_from_vot(xml,linkfile)
    links = read_links(linkfile)
    filenames =  get_files(links,ra,dec,storage,work,False)
    output = work +"output_ra_{}_dec_{}".format(ra,dec)

    make_the_coads(filenames,ra,dec,output,storage)
    os.chdir(work)
    with open("coords_xcs.txt") as f:
        xcs_coords=f.readlines()
    with open("coords_peak.txt") as f:
        peak_coords = f.readlines()
    final_images(xcs_coords,peak_coords,output,455,fileno)
    final_images(xcs_coords,peak_coords,output,909,fileno)
    final_images(xcs_coords,peak_coords,output,1365,fileno)
    final_images(xcs_coords,peak_coords,output,1818,fileno)
