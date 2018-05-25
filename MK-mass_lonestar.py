#module load python/2.7.11

import sys
import emcee
import numpy as np
import time
from time import gmtime, strftime
import scipy.optimize as op

#import pyfits

## nvar fehon smallstep largestep threads

smallstep = 10000
bigstep =  750000
nvar = 6
fehon = 0
erron = 1

if len(sys.argv) >= 2:
    nvar = int(sys.argv[1])
    fehon = int(sys.argv[2])
if len(sys.argv) >= 5:
    smallstep = float(sys.argv[3])
    bigstep = float(sys.argv[4])
threads = 72
if len(sys.argv) >= 6:
    threads = int(sys.argv[5])

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
num_lines = sum(1 for line in open('data2.txt'))
nbin = num_lines-1 ## number of binaries in sample
print nvar,fehon,erron,smallstep,bigstep,threads
print 'Number of targets:',nbin

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
mcsize = 10000
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


#result_ben = np.array([0.2311,-0.1352, 0.0400, 0.0038, -0.0032]) # benedict fit value
#result1 = np.array([0.23323026,-0.10887911, 0.019990399, 0.00027286744, -0.00046073982])# Mann fit value
#result_delf = [0.001*1.8,0.001*6.12,0.001*13.205,-6.2315*0.001,0.001*0.37529]
result2 = plxval
errfac = np.array([np.log(0.01)])
result3_base = np.array([ -0.647825172838, -0.209752211589,   0.00198928002858,  0.00559357629321, 0.000109819339823])
result3_base = np.array([-0.64494970,-0.20761498,-0.0038724139,0.0055028137,0.00047754176])
result3_base = np.array([ -0.64871221,-0.21016848,0.0023172133,0.0056601840,-0.00012339957])
fehcoeff = np.array([0.01])
if nvar-fehon-erron == 4:
    result3 = np.array([-0.64929964,-0.19901594,0.0011934525,0.0032461656])
if nvar-fehon-erron == 5:
    result3 = np.array([ -0.64871221,-0.21016848,0.0023172133,0.0056601840,-0.00012339957])
if nvar-fehon-erron == 6:
    result3 = np.array([-0.64811677,-0.20157744, 0.00017210980,0.0040999402,-0.00011201641,1e-6])
if nvar-fehon-erron == 7:
    result3 = np.array([-0.64811677,-0.20157744, 0.00017210980,0.0040999402,-0.00011201641,1e-6,1e-6])
if erron == 1:
    if fehon == 1:
        result = np.concatenate([result3,fehcoeff,errfac,result2])
    if fehon != 1:
        result = np.concatenate([result3,errfac,result2])
if erron != 1:
    if fehon == 1:
        result = np.concatenate([result3,fehcoeff,result2])
    if fehon != 1:
        result = np.concatenate([result3,result2])
    
# In[ ]:

factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
empmass = factor*smaper/plxval**3
e_empmass = empmass*np.sqrt((esmaper/smaper)**2 +9.0*(plxprior/plxval)**2)


# In[ ]:

## this is mostly for checking things are reasonable
mka = kp - 5.0*(np.log10(1000.0/plxval)-1.)
mkb = ks - 5.0*(np.log10(1000.0/plxval)-1.)
a, b, c, d, e = result3_base
mka_err = ekp
mkb_err = eks
mass1 = 10.0**(a + b*(mka-7.5) + c*(mka-7.5)**2 + d*(mka-7.5)**3 + e*(mka-7.5)**4)
mass2 = 10.0**(a + b*(mkb-7.5) + c*(mkb-7.5)**2 + d*(mkb-7.5)**3 + e*(mkb-7.5)**4)
mass1_err = (np.log(10)*(b+2*c*(mka-7.5)+3*d*(mka-7.5)**2+4*e*(mka-7.5)**3))*mass1*mka_err
mass2_err = (np.log(10)*(b+2*c*(mkb-7.5)+3*d*(mkb-7.5)**2+4*e*(mkb-7.5)**3))*mass2*mkb_err

model_err = np.sqrt(mass1_err**2+mass2_err**2)
model_err2 = np.sqrt(mass1_err**2+mass2_err**2)#+(scat*mass1)**2+(scat*mass2)**2 )
model = mass1+mass2

for i in range(0,len(empmass)):
    print "{:16s}".format(name[i]),     "{0:.3f}".format(empmass[i]),"{0:.3f}".format(e_empmass[i]),     "{0:.4f}".format(model[i]),"{0:.4f}".format(model_err[i]),     "{0:.4f}".format(mka[i]),"{0:.4f}".format(mkb[i]),     "{0:.3f}".format(mass1[i]),"{0:.3f}".format(mass2[i]),     "{0:.1f}".format(np.abs(empmass[i]-model[i])/np.sqrt(e_empmass[i]**2+model_err[i]**2)),     "{0:.2f}".format(feh[i])#     "{0:.3f}".format(ekp[i]),"{0:.3f}".format(eks[i]),
    
print 'rough rchi^2:',np.sum((empmass-model)**2/(e_empmass**2+model_err**2))/(empmass.size-5.)


# In[ ]:

def lnlike(theta, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon, erron):
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
    for ii in range(nvar-fehon-erron):
        factor1 += theta[ii]*mka**(ii)
        factor2 += theta[ii]*mkb**(ii)
    mass1 = 10.0**(factor1)
    mass2 = 10.0**(factor2)
    if fehon == 1:
        mass1*=(1+theta[nvar-1-erron]*feh)
        mass2*=(1+theta[nvar-1-erron]*feh)
    if np.min(mass1) <= 0 or np.min(mass2) <= 0:
        return -np.inf
    
    ## this is where we check to see if the relation always does brighter=higher mass (if not return -np.inf)
    mk = np.linspace(4.0,11.1,100) - zp
    tmp = mk*0
    for ii in range(nvar-fehon-erron):
        tmp += theta[ii]*mk**ii
    l = 10.0**tmp
    check = all(l[i] >= l[i+1] for i in xrange(len(l)-1))
    if not check:
        return -np.inf

    mka_err = ekp
    mkb_err = eks
    factor1 = mka_err*0
    factor2 = mkb_err*0
    for ii in range(nvar-fehon-erron):
      factor1 += ii*theta[ii]*mka**(ii-1)
      factor2 += ii*theta[ii]*mkb**(ii-1)
    mass1_err = np.abs((np.log(10.)*(factor1))*mass1*mka_err)
    mass2_err = np.abs((np.log(10.)*(factor2))*mass2*mkb_err)
    
    model_err = np.sqrt(mass1_err**2+mass2_err**2)
    model = mass1+mass2
    if erron == 1:
        lnf = theta[nvar-1]
        inv_sigma2 = 1.0/(e_empmass**2+model_err**2 + model**2*np.exp(2*lnf))
    if erron != 1:
        inv_sigma2 = 1.0/(e_empmass**2+model_err**2)
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

def lnprob(theta, plxval, plxprior, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon, erron):
    lp = lnprior(theta, plxval, plxprior, nvar)
    if not np.isfinite(lp):
        return -np.inf
    like = lnlike(theta, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon, erron)
    if not np.isfinite(like):
        return -np.inf
    val = lp + like
    return val


# In[ ]:

## some tests:
theta = result
print 'make sure this is not -inf:'
lnl = lnlike(theta, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon, erron)
print lnl
if ~np.isfinite(lnl):
    sys.exit("bad starting conditions")

start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", gmtime())
ndim, nwalkers = result.size, 400
pos = [result + 1e-2*result*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, 
                                args=(plxval, plxprior, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon, erron),
                               threads=threads)
## burn-in and/or testing
pos, prob, state = sampler.run_mcmc(pos, smallstep)
print 'Finished burn/test phase runtime (minutes):'
print (time.time() - start_time)/60
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
print 'Expected time (hours):'
print ((time.time() - start_time)/60)*(bigstep/smallstep)/60
print time.strftime("%a, %d %b %Y %H:%M:%S", gmtime())

dat = sampler.chain
prob = sampler.lnprobability
accept = sampler.acceptance_fraction
dat = sampler.flatchain
like = sampler.flatlnprobability
arr = dat
best = (like == max(like))
a = arr[best,0:nvar]
print a[0,:]
print np.mean(arr[:,nvar-1]), np.std(arr[:,nvar-1])

sampler.reset()
start_time = time.time()
nsteps = bigstep
thin = 500
kwargs = {'thin': thin}
print 'Starting run!'
sampler.sample(pos, iterations=nsteps, **kwargs)
#for i, result in enumerate(sampler.sample(pos, iterations=nsteps, **kwargs)):
#    if (i+1) % 50000 == 0:
#        print("{0:5.1%}".format(float(i) / nsteps)),
#        ("{0:5.2%}".format((time.time() - start_time)/60))
print ' '
print 'Done, runtime (hours):'
print (time.time() - start_time)/3600
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

short = sampler.chain[:,:,0:nvar]
adder = ''
if fehon == 1: adder = '_feh'
## save out the relevant chains
if bigstep > 10000:
    np.save('/work/03344/amann/lonestar/ML/Mk-M_'+str(nvar)+adder+'_short', np.array(short))
    np.save('/work/03344/amann/lonestar/ML/Mk-M_'+str(nvar)+adder+'_chain', np.array(sampler.chain))
    np.save('/work/03344/amann/lonestar/ML/Mk-M_'+str(nvar)+adder+'_flat', np.array(sampler.flatchain))
    np.save('/work/03344/amann/lonestar/ML/Mk-M_'+str(nvar)+adder+'_lnprob', np.array(sampler.lnprobability))
    np.save('/work/03344/amann/lonestar/ML/Mk-M_'+str(nvar)+adder+'_accept',np.array(sampler.acceptance_fraction))

#pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'_short.fits', short, clobber=True)
#pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'.fits', sampler.chain, clobber=True)
#pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'_accept.fits', sampler.acceptance_fraction, clobber=True)
#pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'_lnprob.fits', sampler.lnprobability, clobber=True)


