import numpy, aplpy
import astropy.io.fits as pyfits
 
objno=numpy.loadtxt('temp_objno.txt',dtype=str)
objno=str(objno)
tpt,objno=objno.split('/')
ra,dec=objno.split('_')
cdlist=numpy.loadtxt('forcropandcontours.txt',dtype=str)
# must be index,ra,dec,obsid

for i in range(len(cdlist)):
  if str(cdlist[i][1])==ra:
    racent=cdlist[i][1]
    deccent=cdlist[i][2]
    xrayname=cdlist[i][3]
print ra,dec,xrayname
# CREATE 3-BAND COLOUR IMAGE

#g_data=pyfits.open('tempfits/coadd_g.fits')[0].data
#r_data=pyfits.open('tempfits/coadd_r.fits')[0].data
#i_data=pyfits.open('tempfits/coadd_i.fits')[0].data
#vmg=numpy.median(g_data)+numpy.std(g_data)
#vmr=numpy.median(r_data)+numpy.std(r_data)
#vmi=numpy.median(i_data)+numpy.std(i_data)
#vmig=numpy.median(g_data)
#vmir=numpy.median(r_data)
#vmii=numpy.median(i_data)

aplpy.make_rgb_cube(['tempfits/coadd_i.fits','tempfits/coadd_r.fits','tempfits/coadd_g.fits'],'tempfits/opt.fits')   
oi=aplpy.FITSFigure('tempfits/opt.fits',dimensions=[0,1],slices=[0])
#aplpy.make_rgb_image('tempfits/opt.fits','tempfits/opt2.png',stretch_r='arcsinh',stretch_g='arcsinh',stretch_b='arcsinh',vmax_r=vmg,vmax_g=vmr,vmax_b=vmi,vmin_r=vmig,vmin_g=vmir,vmin_b=vmii)
aplpy.make_rgb_image('tempfits/opt.fits','tempfits/opt.png',stretch_r='sqrt',stretch_g='sqrt',stretch_b='sqrt')

#OVERLAY XRAY CONTOURS AND SAVE
try:
  if len(xrayname)==9:
    oi.show_contour('fits_xcs/0%s-0.50-2.00keVmerged_img.fits'%xrayname,overlap=True,levels=[0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6],colors='white',smooth=3)
  elif len(xrayname)==8:
    oi.show_contour('fits_xcs/00%s-0.50-2.00keVmerged_img.fits'%xrayname,overlap=True,levels=[0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6],colors='white',smooth=3)
  elif len(xrayname)==7:
    oi.show_contour('fits_xcs/000%s-0.50-2.00keVmerged_img.fits'%xrayname,overlap=True,levels=[0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6],colors='white',smooth=3)
except: pass

#oi.show_rgb('images/%s_%s_c.jpeg'%(ra,dec)) #use SDSS image instead
oi.show_rgb('tempfits/opt.png')
oi.axis_labels.hide()
oi.tick_labels.hide()
oi.ticks.hide()
oi.set_theme('publication')

imgsize=12.
oi.recenter(float(racent),float(deccent),radius=imgsize/(2*60.))
oi.save('images/%s_%s_c.png'%(racent,deccent),adjust_bbox=True,dpi=100)

imgsize=6.
oi.recenter(float(racent),float(deccent),radius=imgsize/(2*60.))
oi.save('images/%s_%s_b.png'%(racent,deccent),adjust_bbox=True,dpi=100)

imgsize=3.
oi.recenter(float(racent),float(deccent),radius=imgsize/(2*60.))
oi.save('images/%s_%s_a.png'%(racent,deccent),adjust_bbox=True,dpi=100)

oi.close() 
