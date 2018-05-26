# program converting the .npy into .fits files

import sys
import numpy as np
import pyfits
import glob
import matplotlib
import corner
from matplotlib.backends.backend_pdf import PdfPages

list = (glob.glob("/Users/andrewmann/Dropbox/Post_ML/*.npy"))

for i in list:
    tmp = np.load(i)
    file = i.replace('.npy','.fits')
    print file
    pyfits.writeto(file, tmp, clobber=True)


## this part makes plots
#list = (glob.glob("/Users/andrewmann/Dropbox/Post_ML/*flat.npy"))
list = (glob.glob("/Users/andrewmann/Dropbox/Post_ML/*flat.npy"))
#list2 = (glob.glob("/Users/andrewmann/Dropbox/Post_ML/*s.fits.npy"))
#list = [list1,list2]
#print list
for i in list:
    flat = np.load(i)
    print i
    pdfname = 'fail'
    if i.find("feh") != -1:
        if i.find("8") != -1:
            labels = [r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$a_5$',r'$f$',r'$e$']
            pdfname = 'output_params_8feh.pdf'
            flat = flat[:,0:8]
            flat[:,7] = np.exp(flat[:,7])/2.0
        if i.find("7") != -1:
            labels = [r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$f$',r'$e$']
            pdfname = 'output_params_7feh.pdf'
            flat = flat[:,0:7]
            flat[:,6] = np.exp(flat[:,6])/2.0
        if i.find("6") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$f$',r'$e$']
            pdfname = 'output_params_6feh.pdf'
            flat = flat[:,0:6]
            flat[:,5] = np.exp(flat[:,5])/2.0
    else:
        if i.find("7") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$a_5$',r'$e$']
            pdfname = 'output_params_7.pdf'
            flat = flat[:,0:7]
            flat[:,6] = np.exp(flat[:,6])/2.0
        if i.find("6") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$e$']
            pdfname = 'output_params_6.pdf'
            flat = flat[:,0:6]
            flat[:,5] = np.exp(flat[:,5])/2.0
        if i.find("5") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$e$']
            pdfname= 'output_params_5.pdf'
            flat = flat[:,0:5]
            flat[:,4] = np.exp(flat[:,4])/2.0
    print flat.shape
    print labels,pdfname
    fig = corner.corner(flat, labels=labels, show_titles=True, title_kwargs={"fontsize": 11},title_fmt='.4f',quantiles=(0.16, 0.84), levels=[(1-np.exp(-0.5)),(1-np.exp(-2)),(1-np.exp(-4.5))])
    pp = PdfPages(pdfname)
    pp.savefig(fig)
    pp.close()
