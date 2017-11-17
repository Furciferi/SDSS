import sys, re, os, aplpy, numpy
import astropy.io.fits as pyfits
counter = 0

table = pyfits.open('')
data = table[1].data

NAME = data.field('NAME')
RA = data.field('RA_XCS')
DEC = data.field('DEC_XCS')
ID = data.field('OBSID')
THUMBNAME = data.field('THUMBNAME')
RA_2 = data.field('RA_PEAK')
DEC_2 = data.field('DEC_PEAK')
RA_3 = data.field('RA_RM')
DEC_3 = data.field('DEC_RM')

coordAxes = False

# Making STIFF Images
#for counter in range(int(sys.argv[1]),int(sys.argv[2])):
for counter in range(len(NAME)):
  print('in = '+str(counter))

  oi = aplpy.FITSFigure(THUMBNAME[counter]+'_i.fits',dimensions=[0,1],slices=[0])
  oi.show_rgb(str(THUMBNAME[counter]+'.jpg'))
  oi.axis_labels.hide()
  oi.tick_labels.hide()
  oi.ticks.hide()

# Adding markers to the positions on the image
  oi.show_markers(RA, DEC, edgecolor='cyan',linestyle='--',marker="o",s=500)
  oi.show_markers(RA_2, DEC_2, edgecolor='magenta',linestyle='--',marker="o",s=500)
  oi.show_markers(RA_3, DEC_3, edgecolor='yellow',linestyle='--',marker="o",s=500)

  print('saving...')

# if the size of the total cutout  is 9' x 9'
# 1' x 1'
  oi.recenter(float(RA[counter]),float(DEC[counter]),radius=0.01666)
  oi.save(str(NAME[counter])+'_'+str(ID[counter])+'_ctr_d.png',dpi=100)
# 3' x 3'
  oi.recenter(float(RA[counter]),float(DEC[counter]),radius=0.025)
  oi.save(str(NAME[counter])+'_'+str(ID[counter])+'_ctr_a.png',dpi=100)
# 6' x 6'
  oi.recenter(float(RA[counter]),float(DEC[counter]),radius=0.05)
  oi.save(str(NAME[counter])+'_'+str(ID[counter])+'_ctr_b.png',dpi=100)
# 9' x 9'
  oi.recenter(float(RA[counter]),float(DEC[counter]),radius=0.075)
  oi.save(str(NAME[counter])+'_'+str(ID[counter])+'_ctr_c.png',dpi=100)

  oi.close()
