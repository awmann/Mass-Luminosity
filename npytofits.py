# program converting the .npy into .fits files

import sys
import numpy as np
import pyfits
import glob
import matplotlib
import corner
import os.path
from matplotlib.backends.backend_pdf import PdfPages

list = (glob.glob("../Post_ML/*.npy"))
#list = (glob.glob("../Post_ML/*feh2*.npy"))

for i in list:
    file = i.replace('.npy','.fits')
    if not os.path.isfile(file):
        print file
        tmp = np.load(i)
        pyfits.writeto(file, tmp, clobber=True)


## this part makes plots
list = (glob.glob("../Post_ML/*flat.npy"))
#list = (glob.glob("../*feh2*er1*flat.npy"))
#list2 = (glob.glob("../Post_ML/*s.fits.npy"))
#list = [list1,list2]
#print list

## how do we burn-in this shit?
for i in list:
    flat = np.load(i)
    print i
    pdfname = 'fail'
    adder = '?'
    dif = 1.0
    if i.find("er1") != -1:
        adder = '_eMs'
        dif = 1.0#np.sqrt(2)
    if i.find("er2") != -1:
        adder = '_eplx'
        dif = np.sqrt(2)/3.
    if i.find("er3") != -1:
        adder = '_eMK'
        dif = 1.0
    if i.find("feh2_") != -1:
        if i.find("10") != -1:
            labels = [r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$a_5$',r'$a_6$',r'$f$',r'$g$',r'$\sigma_e$']
            pdfname = 'output_10feh2'+adder+'.pdf'
            flat = flat[:,0:10]
            flat[:,9] = np.exp(flat[:,9])/dif
        if i.find("9") != -1:
            labels = [r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$a_5$',r'$f$',r'$g$',r'$\sigma_e$']
            pdfname = 'output_9feh2'+adder+'.pdf'
            flat = flat[:,0:9]
            flat[:,8] = np.exp(flat[:,8])/dif
        if i.find("8") != -1:
            labels = [r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$f$',r'$g$',r'$\sigma_e$']
            pdfname = 'output_8feh2'+adder+'.pdf'
            flat = flat[:,0:8]
            flat[:,7] = np.exp(flat[:,7])/dif
        if i.find("7") != -1:
            labels = [r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$f$',r'$g$',r'$\sigma_e$']
            pdfname = 'output_7feh2'+adder+'.pdf'
            flat = flat[:,0:7]
            flat[:,6] = np.exp(flat[:,6])/dif
    if i.find("feh_") != -1:
        if i.find("9") != -1:
            labels = [r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$a_5$',r'$a_6$',r'$f$',r'$\sigma_e$']
            pdfname = 'output_9feh'+adder+'.pdf'
            flat = flat[:,0:9]
            flat[:,8] = np.exp(flat[:,8])/dif
        if i.find("8") != -1:
            labels = [r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$a_5$',r'$f$',r'$\sigma_e$']
            pdfname = 'output_8feh'+adder+'.pdf'
            flat = flat[:,0:8]
            flat[:,7] = np.exp(flat[:,7])/dif
        if i.find("7") != -1:
            labels = [r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$f$',r'$\sigma_e$']
            pdfname = 'output_7feh'+adder+'.pdf'
            flat = flat[:,0:7]
            flat[:,6] = np.exp(flat[:,6])/dif
        if i.find("6") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$f$',r'$\sigma_e$']
            pdfname = 'output_6feh'+adder+'.pdf'
            flat = flat[:,0:6]
            flat[:,5] = np.exp(flat[:,5])/dif
    if i.find("feh") == -1:
        if i.find("8") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$a_5$',r'$a_6$',r'$\sigma_e$']
            pdfname = 'output_8'+adder+'.pdf'
            flat = flat[:,0:8]
            flat[:,7] = np.exp(flat[:,7])/dif
        if i.find("7") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$a_5$',r'$\sigma_e$']
            pdfname = 'output_7'+adder+'.pdf'
            flat = flat[:,0:7]
            flat[:,6] = np.exp(flat[:,6])/dif
        if i.find("6") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$\sigma_e$']
            pdfname = 'output_6'+adder+'.pdf'
            flat = flat[:,0:6]
            flat[:,5] = np.exp(flat[:,5])/dif
        if i.find("5") != -1:
            labels=[r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$',r'$\sigma_e$']
            pdfname= 'output_5'+adder+'.pdf'
            flat = flat[:,0:5]
            flat[:,4] = np.exp(flat[:,4])/dif
    print flat.shape
    print labels,pdfname
    fig = corner.corner(flat, labels=labels, show_titles=True, title_kwargs={"fontsize": 14},title_fmt='.5f',quantiles=(0.16, 0.84), levels=[(1-np.exp(-0.5)),(1-np.exp(-2)),(1-np.exp(-4.5))])
    pp = PdfPages(pdfname)
    pp.savefig(fig)
    pp.close()
