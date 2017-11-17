import numpy
import astropy.io.fits as pyfits
from astropy import wcs 
import pyraf
from pyraf import iraf

objno=numpy.loadtxt('temp_objno.txt',dtype=str)
objno=str(objno)
tpt,objno=objno.split('/')
ra,dec=objno.split('_')
cdlist=numpy.loadtxt('forcropandcontours.txt',dtype=str)
# must be index,ra,dec,obsid
filelist=numpy.loadtxt('temp_filesindir.txt',dtype=str) 

for i in range(len(cdlist)):
  if str(cdlist[i][1])==ra:
    racent=cdlist[i][1]
    deccent=cdlist[i][2]
    xrayname=cdlist[i][3]
print '---------------ra and dec are', racent, ra, deccent, dec
imgsize=12
#CUT FITS TO 12'X12'
radeccent=numpy.array([[racent,deccent]],numpy.float_)
for i in range(len(filelist)):
  hdulist=pyfits.open('%s'%filelist[i])
  centre=wcs.WCS(hdulist[0].header).wcs_world2pix(radeccent,1)
  pxscale=numpy.sqrt((hdulist[0].header['CD1_1'])**2+(hdulist[0].header['CD2_1'])**2) #in degrees
  imgpxsize=(imgsize/60.)/pxscale #12 arcmin in pixels
  if int(centre[0][0]-imgpxsize/2.)<0 : 
    xst=1#0
  else:
    xst=int(centre[0][0]-imgpxsize/2.)
  if int(centre[0][0]+imgpxsize/2.)>int(hdulist[0].header['NAXIS1']) : 
    xen=int(hdulist[0].header['NAXIS1'])
  else:
    xen=int(centre[0][0]+imgpxsize/2.)
  if int(centre[0][1]-imgpxsize/2.)<0 : 
    yst=1#0
  else:
    yst=int(centre[0][1]-imgpxsize/2.)
  if int(centre[0][1]+imgpxsize/2.)>int(hdulist[0].header['NAXIS2']) : 
    yen=int(hdulist[0].header['NAXIS2'])
  else:
    yen=int(centre[0][1]+imgpxsize/2.)
  print '---------------xpixel and ypixel are',centre[0][0],centre[0][1]
  pyfits.writeto('tempfits/temp_%s%s.fits'%(str(filelist[i][15]),i),header=hdulist[0].header,data=hdulist[0].data)
  if (xen>0 and yen>0 and xen>xst and yen>yst):  
    iraf.imutil.imcopy(input='tempfits/temp_%s%s.fits['%(str(filelist[i][15]),i)+str(xst)+':'+str(xen)+','+str(yst)+':'+str(yen)+']',output='tempfits/%s%s.fits'%(str(filelist[i][15]),i))
