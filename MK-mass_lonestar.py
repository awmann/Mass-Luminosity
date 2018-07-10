#module load python/2.7.11

import sys
import emcee
import numpy as np
import time
from time import gmtime, strftime
import scipy.optimize as op
import os

## to do:
## reload option (pick up where left off)

#import pyfits

## nvar fehon smallstep largestep threads

smallstep = 10000
bigstep =  750000
nvar = 7
fehon = 0
erron = 1
walkers = 400

## arg = nvar fehon erron smallstep bigstep threads

if len(sys.argv) >= 2:
    nvar = int(sys.argv[1])
    fehon = int(sys.argv[2])
if len(sys.argv) >= 5:
    erron = int(sys.argv[3])
if len(sys.argv) >= 6:
    smallstep = float(sys.argv[4])
    bigstep = float(sys.argv[5])
threads = 72
if len(sys.argv) >= 7:
    threads = int(sys.argv[6])

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
print 'Nvar:',nvar
if fehon == 0: print 'no [Fe/H] term'
if fehon == 1: print '[Fe/H] included'
if fehon == 2: print '[Fe/H] and [Fe/H]*mk term'
if erron == 0:
    print 'no extra error term'
if erron == 1:
    print 'error added to model mass'
if erron == 2:
    print 'error added to plx'
if erron == 3:
    print 'error added to K mags'
print 'burn-in:',smallstep
print 'steps:',bigstep
print 'cpu threads:',threads
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
    

errsub = 0
if erron > 0:
   errsub = 1

#result_ben = np.array([0.2311,-0.1352, 0.0400, 0.0038, -0.0032]) # benedict fit value
#result1 = np.array([0.23323026,-0.10887911, 0.019990399, 0.00027286744, -0.00046073982])# Mann fit value
#result_delf = [0.001*1.8,0.001*6.12,0.001*13.205,-6.2315*0.001,0.001*0.37529]
result2 = plxval
#result3_base = np.array([ -0.647825172838, -0.209752211589,   0.00198928002858,  0.00559357629321, 0.000109819339823])
#result3_base = np.array([-0.64494970,-0.20761498,-0.0038724139,0.0055028137,0.00047754176])
#result3_base = np.array([ -0.652530,-0.203478,0.004481558,0.005100981,-0.000298238])
#result3_base = np.array([ -0.646609, -0.212461,  -0.002653431, 0.007948519, 0.0003689931,  -0.0001922619])
result3_base = np.array([ -0.64831778, -0.21504067,  -7.0397849e-5, 0.0088852845, 0.00013451022,  -0.00024134927],dtype=np.float128)
if erron == 1:
    errfac = np.array([np.log(0.01)],dtype=np.float128)
if erron == 2:
    errfac = np.array([np.log(0.007)],dtype=np.float128)
if erron == 3:
    errfac = np.array([np.log(0.070)],dtype=np.float128)
    
fehcoeff = np.array([0.037],dtype=np.float128)
fehcoeff2 = np.array([-0.041],dtype=np.float128)
if nvar-fehon-errsub == 4:
    result3 = np.array([-0.646460,-0.196306,0.00007207469,0.004148834],dtype=np.float128)
if nvar-fehon-errsub == 5:
    result3 = np.array([ -0.65253050,-0.20347828,0.0044815579,0.0051009810,-0.00029823842],dtype=np.float128)
if nvar-fehon-errsub == 6:
    result3 = np.array([-0.64835236,-0.21556456, -0.00084600047,0.0083775937,0.00016819974,-0.00018743871],dtype=np.float128)
if nvar-fehon-errsub == 7:
    result3 = np.array([-0.64351569,-0.20943173, -0.0076165544,0.0061251049,0.0015557089,-3.0257191e-5,-7.0698029e-5],dtype=np.float128)
#print result3
if erron >= 1:
    if fehon == 2:
        result = np.concatenate([result3,fehcoeff,fehcoeff2,errfac,result2])
    if fehon == 1:
        result = np.concatenate([result3,fehcoeff,errfac,result2])
    if fehon == 0 :
        result = np.concatenate([result3,errfac,result2])
if erron == 0:
    if fehon == 2:
        result = np.concatenate([result3,fehcoeff,fehcoeff2,result2])
    if fehon == 1:
        result = np.concatenate([result3,fehcoeff,result2])
    if fehon == 0:
        result = np.concatenate([result3,result2])

# In[ ]:

factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
empmass = factor*smaper/plxval**3
e_empmass = empmass*np.sqrt((esmaper/smaper)**2 +9.0*(plxprior/plxval)**2)


# In[ ]:

## this is mostly for checking things are reasonable
mka = kp - 5.0*(np.log10(1000.0/plxval)-1.)
mkb = ks - 5.0*(np.log10(1000.0/plxval)-1.)
a, b, c, d, e, aa = result3_base
mka_err = ekp
mkb_err = eks
mass1 = 10.0**(a + b*(mka-7.5) + c*(mka-7.5)**2. + d*(mka-7.5)**3. + e*(mka-7.5)**4. + aa*(mka-7.5)**5.)
mass2 = 10.0**(a + b*(mkb-7.5) + c*(mkb-7.5)**2. + d*(mkb-7.5)**3. + e*(mkb-7.5)**4. + aa*(mkb-7.5)**5.)
mass1_err = (np.log(10)*(b+2*c*(mka-7.5)+3*d*(mka-7.5)**2+4*e*(mka-7.5)**3+5*aa*(mka-7.5)**4))*mass1*mka_err
mass2_err = (np.log(10)*(b+2*c*(mkb-7.5)+3*d*(mkb-7.5)**2+4*e*(mkb-7.5)**3+5*aa*(mkb-7.5)**4))*mass2*mkb_err
    
model_err = np.sqrt(mass1_err**2+mass2_err**2)
model_err2 = np.sqrt(mass1_err**2+mass2_err**2)
model = mass1+mass2

print 'Name  empmass err  modelmass  err   mka   mkb  mass1 mass2 diff(sig)  feh '
for i in range(0,len(empmass)):
    print "{:16s}".format(name[i]),     "{0:.3f}".format(empmass[i]),"{0:.3f}".format(e_empmass[i]),     "{0:.4f}".format(model[i]),"{0:.4f}".format(model_err[i]),     "{0:.4f}".format(mka[i]),"{0:.4f}".format(mkb[i]),     "{0:.3f}".format(mass1[i]),"{0:.3f}".format(mass2[i]),     "{0:.1f}".format(np.abs(empmass[i]-model[i])/np.sqrt(e_empmass[i]**2+model_err[i]**2)),     "{0:.2f}".format(feh[i])#     "{0:.3f}".format(ekp[i]),"{0:.3f}".format(eks[i])
print 'rough rchi^2:',np.sum((empmass-model)**2/(e_empmass**2+model_err**2))/(empmass.size-5.)

print 'mean diff, weighted mean diff: ',np.mean(empmass-model),np.average(empmass-model,weights=(1./np.sqrt(e_empmass**2+model_err**2)))


def lnlike(theta, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon, erron):
    errsub = 0
    if erron > 0:
        errsub = 1
        lnf = theta[nvar-1]
        if np.exp(lnf) < 0.002:
            return -np.inf

    zp = 7.5
    au = 1.496e13
    msun = 1.989e33
    g = 6.6743e-8 
    mplx = theta[nvar:theta.size]

    if np.min(mplx) <= 0:
        return -np.inf
    factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
    empmass = factor*smaper/mplx**3
    e_empmass = empmass*(esmaper/smaper)
    #if erron == 2: ## error as extra parallax (but in mass only?)
    #    e_empmass = empmass*np.sqrt((esmaper/smaper)**2 + 9.0*(np.exp(theta[nvar-1])**2))
    
    mka = kp - 5.0*(np.log10(1000.0/mplx)-1.) - zp
    mkb = ks - 5.0*(np.log10(1000.0/mplx)-1.) - zp
    if np.min(mka) <= -9:
        #print np.min(mka)
        return -np.inf
    factor1 = mka*0
    factor2 = mkb*0
    for ii in range(nvar-fehon-errsub):
        factor1 += theta[ii]*mka**(ii)
        factor2 += theta[ii]*mkb**(ii)
        
    if np.max(factor1) > 1:
        return -np.inf
    if np.max(factor2) > 1:
        return -np.inf
        
    mass1 = 10.0**(factor1)
    mass2 = 10.0**(factor2)
    mass1_base = mass1
    mass2_base = mass2
    if np.min(mass1_base) < 0 or np.min(mass2_base) < 0 or np.max(mass1_base) > 1.0 or np.max(mass2_base) > 1.0:
        return -np.inf
        
    if fehon == 1:
        f = theta[nvar-fehon-errsub]
        g = 0.0
        mass1*=(1.+f*feh)
        mass2*=(1.+f*feh)
    if fehon == 2:
        f = theta[nvar-fehon-errsub]
        g = theta[nvar-fehon-errsub+1]
        mass1*=(1.+(f+g*mka)*feh)
        mass2*=(1.+(f+g*mkb)*feh)
    if np.min(mass1) <= 0 or np.min(mass2) <= 0:
        return -np.inf

    if np.min(mass1) < 0 or np.min(mass2) < 0 or np.max(mass1) > 1.0 or np.max(mass2) > 1.0:
        return -np.inf
    ## this is where we check to see if the relation always does brighter=higher mass (if not return -np.inf)
    mk = np.linspace(4.0,11.5,100) - zp
    tmp = mk*0
    for ii in range(nvar-fehon-errsub):
        tmp += theta[ii]*mk**ii
    l = 10.0**tmp
    check = all(l[i] >= l[i+1] for i in xrange(len(l)-1))
    if not check:
        return -np.inf

    mka_err = ekp
    mkb_err = eks
    if erron == 3: ## add error to absolute magnitude
        mka_err = np.sqrt(ekp**2+(np.exp(2.*lnf)))
        mkb_err = np.sqrt(eks**2+(np.exp(2.*lnf)))
    if fehon >= 1:
        e_feh = 0.04# global 0.04 (relative!) uncertainties on [Fe/H]
    else:
        e_feh = 0.0
    factor1 = mka_err*0
    factor2 = mkb_err*0
    for ii in range(nvar-fehon-erron):
      factor1 += ii*theta[ii]*mka**(ii-1)
      factor2 += ii*theta[ii]*mkb**(ii-1)
    if fehon == 0:
        mass1_err = np.abs((np.log(10.)*(factor1))*mass1*mka_err)
        mass2_err = np.abs((np.log(10.)*(factor2))*mass2*mkb_err)
    if fehon >= 1:
        mass1_err = np.sqrt((mass1_base*(g*feh+np.log(10.)*factor1*(1.+f*feh+g*mka*feh)))**2*mka_err**2 + ((f+g*mka)*mass1_base)**2*e_feh**2) 
        mass2_err = np.sqrt((mass2_base*(g*feh+np.log(10.)*factor2*(1.+f*feh+g*mkb*feh)))**2*mkb_err**2 + ((f+g*mkb)*mass2_base)**2*e_feh**2)

    model_err = np.sqrt(mass1_err**2+mass2_err**2)
    model = mass1+mass2
    if erron == 1:
        inv_sigma2 = 1.0/(e_empmass**2+model_err**2 + model**2*np.exp(2*lnf))
    if erron != 1:
        inv_sigma2 = 1.0/(e_empmass**2+model_err**2)
    if np.min(inv_sigma2) <= 0:
        return -np.inf
    return -0.5*(np.sum((empmass-model)**2*inv_sigma2 - np.log(inv_sigma2)))

# In[ ]:

def lnprior(theta, plxval, plxprior, nvar, erron):
    mplx = theta[nvar:theta.size]
    lp = 0
    if erron == 2:
        plxprior = np.sqrt(plxprior**2 + (plxval*(np.exp(theta[nvar-1]))**2))
    if np.min(mplx) <= 0:
        return -np.inf
    lp += np.sum(((mplx-plxval)**2)/(plxprior**2))
    lp*=(-0.5)
    if not np.isfinite(lp):
        return 0.0
    return lp


# In[ ]:

def lnprob(theta, plxval, plxprior, smaper, esmaper, kp, ks, ekp, eks, feh, nvar, fehon, erron):
    lp = lnprior(theta, plxval, plxprior, nvar, erron)
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
#if np.isfinite(lnl):
#    sys.exit("bad starting conditions")

start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", gmtime())
ndim, nwalkers = result.size, walkers
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
print np.mean(np.exp(arr[:,nvar-1])), np.std(np.exp(arr[:,nvar-1]))
#print np.mean(np.exp(arr[:,nvar-2])), np.std(np.exp(arr[:,nvar-2]))

sampler.reset()
start_time = time.time()
nsteps = bigstep
thin = 500
kwargs = {'thin': thin}
print 'Starting run!'
for i, result in enumerate(sampler.sample(pos, iterations=nsteps, **kwargs)):
    if (i+1) % 50000 == 0:
        print("{0:5.1%}".format(float(i) / nsteps)),
        ("{0:5.2%}".format((time.time() - start_time)/60))
print ' '
print 'Done, runtime (hours):'
print (time.time() - start_time)/3600
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

short = sampler.chain[:,:,0:nvar]
adder = ''
if fehon == 2:
    adder = '_feh2'
    if erron == 1:
        adder = '_feh2_er1'
    if erron == 2:
        adder = '_feh2_er2'
    if erron == 3:
        adder = '_feh2_er3'
if fehon == 1:
    adder = '_feh'
    if erron == 1:
        adder = '_feh_er1'
    if erron == 2:
        adder = '_feh_er2'
    if erron == 3:
        adder = '_feh_er3'
if fehon == 0:
    adder = ''
    if erron == 1:
        adder = '_er1'
    if erron == 2:
        adder = '_er2'
    if erron == 3:
        adder = '_er3'
## save out the relevant chains
if bigstep > 1000:
    if '/Users/' not in os.getcwd():
        print 'saving in work'
        np.save('/work/03344/amann/lonestar/ML/Mk-M_'+str(nvar)+adder+'_short', np.array(short))
        np.save('/work/03344/amann/lonestar/ML/Mk-M_'+str(nvar)+adder+'_chain', np.array(sampler.chain))
        np.save('/work/03344/amann/lonestar/ML/Mk-M_'+str(nvar)+adder+'_flat', np.array(sampler.flatchain))
        np.save('/work/03344/amann/lonestar/ML/Mk-M_'+str(nvar)+adder+'_lnprob', np.array(sampler.lnprobability))
        np.save('/work/03344/amann/lonestar/ML/Mk-M_'+str(nvar)+adder+'_accept',np.array(sampler.acceptance_fraction))
    else:
        print 'saving in current directory'
        np.save('Mk-M_'+str(nvar)+adder+'_short', np.array(short))
        np.save('Mk-M_'+str(nvar)+adder+'_chain', np.array(sampler.chain))
        np.save('Mk-M_'+str(nvar)+adder+'_flat', np.array(sampler.flatchain))
        np.save('Mk-M_'+str(nvar)+adder+'_lnprob', np.array(sampler.lnprobability))
        np.save('Mk-M_'+str(nvar)+adder+'_accept',np.array(sampler.acceptance_fraction))
#pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'_short.fits', short, clobber=True)
#pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'.fits', sampler.chain, clobber=True)
#pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'_accept.fits', sampler.acceptance_fraction, clobber=True)
#pyfits.writeto('Mk-Mass_log_emcee'+str(nvar-1)+adder+'_lnprob.fits', sampler.lnprobability, clobber=True)
