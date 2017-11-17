dr = "dr13"
ra= 247.14761112
dec= 39.548403
scale=0.396
height =1818
width =1818
link = ("http://skyserver.sdss.org/{}/SkyServerWS/ImgCutout/getjpeg?TaskName=Skyserver.Chart.List&ra={}&dec={}&scale={}&width={}&height={}&opt=").format(dr,ra,dec,scale,width,height)

print link
