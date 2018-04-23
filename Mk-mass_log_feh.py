
# coding: utf-8

# In[ ]:

import sys
import emcee
import matplotlib
import numpy as np
import corner
import scipy.optimize as op

smallstep = 10000
bigstep =  500000

nvar = int(sys.argv[1])
fehon = int(sys.argv[2])
if len(sys.argv) == 4:
    smallstep = float(sys.argv[3])
    bigstep = float(sys.argv[4])

if smallstep >= bigstep:
    tmp = bigstep
    smallstep = bigstep
    bigstep = tmp

#get_ipython().magic(u'matplotlib inline')

# some important variables
g = 6.6743e-8
msun = 1.989e33
au = 1.496e13
pi = 3.14159


# In[ ]:

nbin = 60 ## number of binaries in sample
#nvar = 5
#fehon = 1 ## turn on Fe/H fit
print nvar,fehon,smallstep,bigstep

## read in delK, parallax
delk = np.zeros(nbin)
edelk = np.zeros(nbin)
plxval = np.zeros(nbin)
plxprior = np.zeros(nbin)
name = strs = ['']*nbin
k = np.zeros(nbin)
ek = np.zeros(nbin)
feh = np.zeros(nbin)
f = open('data2.txt','r')
header1 = f.readline()
i=0
for line in f:
    line = line.strip()
    columns = line.split()
    name[i] = columns[0]
    delk[i] = float(columns[1])
    edelk[i] = float(columns[2])
    plxval[i] = float(columns[3])
    plxprior[i] = float(columns[4])
    k[i] = float(columns[5])
    ek[i] = float(columns[6])
    feh[i] = float(columns[8])
    i+=1
f.close()

## now for the sma**3/per**2
f = open('fits.txt','r')
header1 = f.readline()
i=0
smaper = np.zeros(nbin)
esmaper = np.zeros(nbin)
for line in f:
    line = line.strip()
    columns = line.split()
    smaper[i] = float(columns[0])
    esmaper[i] = float(columns[1])
    i+=1
f.close()

fluxratio = 10.0**(delk/2.5)
del_eps = 2.5*np.log10(1.0+1.0/fluxratio)
kp = del_eps+k
ks = kp + delk

# compute MC errors on Kp, Ks
mcsize = 50000
ekp = kp*0.
eks = ks*0.
for i in range(0,len(ks)):
    ktmp = k[i]+ek[i]*np.random.standard_normal(mcsize)
    deltmp = delk[i]+edelk[i]*np.random.standard_normal(mcsize)
    fluxratio = 10.0**(deltmp/2.5)
    del_eps = 2.5*np.log10(1.0+1.0/fluxratio)
    kpt = del_eps+ktmp
    kst = kp[i] + ktmp
    ekp[i] = np.std(kpt)
    eks[i] = np.std(kst)


result_ben = np.array([0.2311,-0.1352, 0.0400, 0.0038, -0.0032]) # benedict fit value
result1 = np.array([0.23323026,-0.10887911, 0.019990399, 0.00027286744, -0.00046073982])# Mann fit value
result2 = plxval
result_delf = [0.001*1.8,0.001*6.12,0.001*13.205,-6.2315*0.001,0.001*0.37529]

result3= np.array([-0.663541,-0.212175 , 0.00594868 ,  0.00641549, -0.000416446]) ## one I measured in IDL
fehcoeff = np.array([0.01])
result_ben_test = np.concatenate([result_ben,np.array([0.00]),result2])
result_nofeh = np.concatenate([result3,result2])
result3_base = result3
if nvar-fehon == 4:
    result3 = np.array([-0.64929964,-0.19901594,0.0011934525,0.0032461656])
if nvar-fehon == 5:
    result3= np.array([-0.64372527,-0.20007154,-0.0037546659,0.0033291741,0.00027382806])
if fehon == 1:
    result = np.concatenate([result3,fehcoeff,result2])
if fehon != 1:
    result = np.concatenate([result3,result2])
print result[0:nvar]
print result[nvar:result.size]

# In[ ]:

factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
empmass = factor*smaper/plxval**3
e_empmass = empmass*np.sqrt((esmaper/smaper)**2 +9.0*(plxprior/plxval)**2)


# In[ ]:

## this is mostly for checking things are reasonable
mka = kp - 5.0*(np.log10(1000.0/plxval)-1.)
mkb = ks - 5.0*(np.log10(1000.0/plxval)-1.)
a, b, c, d, e = result3_base
mka_err = 0.018## for now
mkb_err = 0.018## for now
mass1 = 10.0**(a + b*(mka-7.5) + c*(mka-7.5)**2 + d*(mka-7.5)**3 + e*(mka-7.5)**4)
mass2 = 10.0**(a + b*(mkb-7.5) + c*(mkb-7.5)**2 + d*(mkb-7.5)**3 + e*(mkb-7.5)**4)
mass1_err = (np.log(10)*(b+2*c*(mka-7.5)+3*d*(mka-7.5)**2+4*e*(mka-7.5)**3))*mass1*mka_err
mass2_err = (np.log(10)*(b+2*c*(mkb-7.5)+3*d*(mkb-7.5)**2+4*e*(mkb-7.5)**3))*mass2*mkb_err

model_err = np.sqrt(mass1_err**2+mass2_err**2)
model_err2 = np.sqrt(mass1_err**2+mass2_err**2)#+(scat*mass1)**2+(scat*mass2)**2 )
model = mass1+mass2

for i in range(0,len(empmass)):
    print "{:10s}".format(name[i]),     "{0:.3f}".format(empmass[i]),"{0:.3f}".format(e_empmass[i]),     "{0:.4f}".format(model[i]),"{0:.4f}".format(model_err[i]),"{0:.3f}".format(100*model_err[i]/model[i]),     "{0:.4f}".format(mka[i]),"{0:.4f}".format(mkb[i]),     "{0:.3f}".format(ekp[i]),"{0:.3f}".format(eks[i]),     "{0:.3f}".format(mass1[i]),"{0:.3f}".format(mass2[i]),     "{0:.1f}".format(np.abs(empmass[i]-model[i])/np.sqrt(e_empmass[i]**2+model_err[i]**2)),     "{0:.2f}".format(feh[i])
    
print np.sum((empmass-model)**2/(e_empmass**2+model_err**2))/(empmass.size-5.)


# In[ ]:

def lnlike(theta, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon):
    zp = 7.5e0
    au = 1.496e13
    msun = 1.989e33
    g = 6.6743e-8 
    #a, b, c, d, e, f = theta[0:nvar]
    mplx = theta[nvar:theta.size]
    if np.min(mplx) <= 0:
        return -np.inf
    factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
    empmass = factor*smaper/mplx**3
    e_empmass = empmass*(esmaper/smaper)**2
    mka = kp - 5.0*(np.log10(1000.0/mplx)-1.) - zp
    mkb = ks - 5.0*(np.log10(1000.0/mplx)-1.) - zp
    factor1 = mka*0
    factor2 = mkb*0
    for ii in range(nvar-fehon):
        factor1 += theta[ii]*mka**ii
        factor2 += theta[ii]*mkb**ii
    mass1 = 10.0**(factor1)
    mass2 = 10.0**(factor2)
    if fehon == 1:
        mass1*=(1+theta[nvar-1])
        mass2*=(1+theta[nvar-1])
    #mass1 = (10.0**(a + b*mka + c*mka**2. + d*mka**3. + e*mka**4.))*(1+f*feh)
    #mass2 = (10.0**(a + b*mkb + c*mkb**2. + d*mkb**3. + e*mkb**4.))*(1+f*feh)
    if np.min(mass1) <= 0 or np.min(mass2) <= 0:
        return -np.inf
    
    ## this is where we check to see if the relation always does brighter=higher mass (if not return -np.inf)
    mk = np.linspace(4.2,11.2,100) - zp
    tmp = mk*0
    for ii in range(nvar-fehon):
        tmp += theta[ii]*mk**ii
    l = 10.0**tmp
    #l = 10.0**(a + b*mk + c*mk**2. + d*mk**3. + e*mk**4.)
    check = all(l[i] >= l[i+1] for i in xrange(len(l)-1))
    if not check:
        return -np.inf
    ## there's a better way to do this... probably to check to see if the slope goes negative

    mka_err = ekp
    mkb_err = eks
    factor1 = mka_err*0
    factor2 = mkb_err*0
    for ii in range(nvar-fehon):
        factor1 += ii*theta[ii]*mka**(ii-1)
        factor2 += ii*theta[ii]*mkb**(ii-1)
    mass1_err = np.abs((np.log(10.)*(factor1))*mass1*mka_err)
    mass2_err = np.abs((np.log(10.)*(factor2))*mass2*mkb_err)
    
    model_err = np.sqrt(mass1_err**2+mass2_err**2)
    model = mass1+mass2
    inv_sigma2 = 1.0/np.sqrt(e_empmass**2+model_err**2)
    return -0.5*(np.sum((empmass-model)**2*inv_sigma2 - np.log(inv_sigma2)))


# In[ ]:

def lnprior(theta, plxval, plxprior, nvar):
    mplx = theta[nvar:theta.size]
    lp = 0
    if np.min(mplx) <= 0:
        return -np.inf
    lp += np.sum(((mplx-plxval)**2)/(plxprior)**2)
    lp*=(-0.5)
    if not np.isfinite(lp):
        return 0.0
    return lp


# In[ ]:

def lnprob(theta, plxval, plxprior, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon):
    lp = lnprior(theta, plxval, plxprior, nvar)
    if not np.isfinite(lp):
        return -np.inf
    like = lnlike(theta, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon)
    if not np.isfinite(like):
        return -np.inf
    val = lp + like
    return val


# In[ ]:

## some tests:
theta = result
print 'make sure this is not -inf:'
print lnlike(theta, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon)
#theta = result_nofeh
#print lnlike(theta, smaper, esmaper, kp, ks, ekp, eks, feh, nvar-1, 0)
#theta = result_ben_test
#print lnlike(theta, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon)


# In[ ]:

import time
from time import gmtime, strftime
start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", gmtime())
ndim, nwalkers = result.size, 600
pos = [result + 1e-2*result*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, 
                                args=(plxval, plxprior, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon),
                               threads=6)
## burn-in and/or testing
pos, prob, state = sampler.run_mcmc(pos, smallstep)
print 'Finished burn/test phase runtime (minutes):'
print (time.time() - start_time)/60
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
print 'Expected time (hours):'
print ((time.time() - start_time)/60)*(bigstep/smallstep)/60
print time.strftime("%a, %d %b %Y %H:%M:%S", gmtime())


# In[ ]:

dat = sampler.chain
prob = sampler.lnprobability
accept = sampler.acceptance_fraction
dat = sampler.flatchain#chain
like = sampler.flatlnprobability
arr = dat
best = (like == max(like))
a = arr[best,0:nvar]#np.median(arr[:,0])
print a[0,:]
print np.mean(arr[:,nvar-1]), np.std(arr[:,nvar-1])


# In[ ]:

sampler.reset()
start_time = time.time()
nsteps = bigstep
thin = 100
kwargs = {'thin': thin }
print 'Starting run!'
for i, result in enumerate(sampler.sample(pos, iterations=nsteps, **kwargs)):
    if (i+1) % 50000 == 0:
        print("{0:5.1%}".format(float(i) / nsteps)),
        ("{0:5.2%}".format((time.time() - start_time)/60))
print ' '
print 'Done, runtime (hours):'
print (time.time() - start_time)/3600
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))


# In[ ]:

print sampler.chain.shape
print sampler.lnprobability.shape
short = sampler.chain[:,:,0:nvar]
print short.shape

adder = ''
if fehon == 1: adder = '_feh'
## save out the relevant chains
import pyfits
pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'_short.fits', short, clobber=True)
pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'.fits', sampler.chain, clobber=True)
pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'_accept.fits', sampler.acceptance_fraction, clobber=True)
pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'_lnprob.fits', sampler.lnprobability, clobber=True)
short = sampler
#pyfits.writeto('Mk-Mass_log_emcee_acor2.fits', sampler.acor, clobber=True)


# In[ ]:

#print dat.shape,nwalkers,nsteps,thin
import corner
from matplotlib.backends.backend_pdf import PdfPages
dat = sampler.flatchain#chain
like = sampler.flatlnprobability
arr = dat
if fehon == 1:
    if nvar == 6:
        labels = [r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$a_5$',r'$f$']
        pdfname = 'output_params_5feh.pdf'
    if nvar == 5:
        labels=[r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$f$']
        pdfname = 'output_params_4feh.pdf'
if fehon == 0:
    if nvar == 5:
        labels=[r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$',r'$a_5$']
        pdfname= 'output_params_5.pdf'
    if nvar == 4:
        labels=[r'$a_1$',r'$a_2$',r'$a_3$',r'$a_4$']
        pdfname = 'output_params_4.pdf'

fig = corner.corner(arr[:,0:nvar], labels=labels, show_titles=True, title_kwargs={"fontsize": 11},title_fmt='.4f',quantiles=(0.16, 0.84), levels=[(1-np.exp(-0.5)),(1-np.exp(-2)),(1-np.exp(-4.5))])
pp = PdfPages(pdfname)
pp.savefig(fig)
pp.close()


# In[ ]:

dat = sampler.chain
prob = sampler.lnprobability
accept = sampler.acceptance_fraction
arr = dat.reshape((dat.shape)[0]*(dat.shape)[1],dat.shape[2])
print 'name plx_mcmc plx_err_mcmc plxinput plx_prior diff_sig1 diff_sig2'
for i in range(nvar,dat.shape[2]):
    print "{:10s}".format(name[i-nvar]),     i,    "{0:.4f}".format(np.median(arr[:,i])),"{0:.4f}".format(np.std(arr[:,i])),"{0:.4f}".format(plxprior[i-nvar]),     "{0:.4f}".format(plxval[i-nvar]),"{0:.4f}".format((plxval[i-nvar]-np.median(arr[:,i]))/plxprior[i-nvar]),     "{0:.4f}".format((plxval[i-nvar]-np.median(arr[:,i]))/np.sqrt(plxprior[i-nvar]**2+np.std(arr[:,i])**2))
    


# In[ ]:

factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
mass = factor*smaper/plxval**3

mka = kp - 5.0*(np.log10(1000.0/plxval)-1.)
mkb = ks - 5.0*(np.log10(1000.0/plxval)-1.)

best = (like == max(like))
a = arr[best,0]#np.median(arr[:,0])
b = arr[best,1]#np.median(arr[:,1])
c = arr[best,2]#np.median(arr[:,2])
d = arr[best,3]#np.median(arr[:,3])
e = arr[best,4]#np.median(arr[:,4])
f = arr[best,5]
print a,b,c,d,e,f

mass1 = (10.0**(a + b*(mka-7.5) + c*(mka-7.5)**2 + d*(mka-7.5)**3 + e*(mka-7.5)**4))*(1+f*feh)
mass2 = (10.0**(a + b*(mkb-7.5) + c*(mkb-7.5)**2 + d*(mkb-7.5)**3 + e*(mkb-7.5)**4))*(1+f*feh)
sysmass = mass1+mass2
a, b, c, d, e = result1
mass1 = a + b*(mka-7.5) + c*(mka-7.5)**2 + d*(mka-7.5)**3 + e*(mka-7.5)**4
mass2 = a + b*(mkb-7.5) + c*(mkb-7.5)**2 + d*(mkb-7.5)**3 + e*(mkb-7.5)**4
sysmass_mann = mass1+mass2
a, b, c, d, e = result_ben
mass1 = a + b*(mka-7.5) + c*(mka-7.5)**2 + d*(mka-7.5)**3 + e*(mka-7.5)**4
mass2 = a + b*(mkb-7.5) + c*(mkb-7.5)**2 + d*(mkb-7.5)**3 + e*(mkb-7.5)**4
sysmass_ben = mass1+mass2


# In[ ]:



