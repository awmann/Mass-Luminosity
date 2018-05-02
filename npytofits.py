# program converting the .npy into .fits files

import sys
import numpy as np
import pyfits
import glob
import matplotlib
import corner
from matplotlib.backends.backend_pdf import PdfPages

list = (glob.glob("/Users/andrewmann/Dropbox/Post_ML2/*.npy"))

for i in list:
    tmp = np.load(i)
    file = i.replace('.npy','.fits')
    print file
    pyfits.writeto(file, tmp, clobber=True)


## this part makes plots
#list = (glob.glob("/Users/andrewmann/Dropbox/Post_ML/*flat.npy"))
list = (glob.glob("/Users/andrewmann/Dropbox/Post_ML2/*flat.npy"))
#list2 = (glob.glob("/Users/andrewmann/Dropbox/Post_ML/*s.fits.npy"))
#list = [list1,list2]
#print list
for i in list:
    flat = np.load(i)
    print i
    pdfname = 'fail'
    if i.find("feh") != -1:
        if i.find("6") != -1:
            labels = [r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$f$']
            pdfname = 'output_params_5feh.pdf'
            flat = flat[:,0:6]
            #flat = np.reshape(sample,(3000000,6))
        if i.find("5") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$f$']
            pdfname = 'output_params_4feh.pdf'
            flat = flat[:,0:5]
            #flat = np.reshape(sample,(3000000,5))
    else:
        if i.find("4") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$']
            pdfname= 'output_params_4.pdf'
            flat = flat[:,0:4]
            #flat = np.reshape(sample,(3000000,5))
        if i.find("5") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$']
            pdfname = 'output_params_5.pdf'
            flat = flat[:,0:5]
            #flat = np.reshape(sample,(3000000,5))
    print flat.shape
    print labels,pdfname
    fig = corner.corner(flat, labels=labels, show_titles=True, title_kwargs={"fontsize": 11},title_fmt='.4f',quantiles=(0.16, 0.84), levels=[(1-np.exp(-0.5)),(1-np.exp(-2)),(1-np.exp(-4.5))])
    pp = PdfPages(pdfname)
    pp.savefig(fig)
    pp.close()
