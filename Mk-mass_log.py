
# coding: utf-8

# In[ ]:


import emcee
import matplotlib
import numpy as np
import corner
import scipy.optimize as op

#%matplotlib inline  

# some important variables
g = 6.6743e-8
msun = 1.989e33
au = 1.496e13
pi = 3.14159


# In[ ]:


nbin = 50 ## number of binaries in sample

## read in delK, parallax
delk = np.zeros(nbin)
edelk = np.zeros(nbin)
plxval = np.zeros(nbin)
plxprior = np.zeros(nbin)
name = strs = ['']*nbin
k = np.zeros(nbin)
ek = np.zeros(nbin)
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
result3= [-0.63373414,-0.20122973,-0.0067369013,0.0024262687,-0.00013746096] ## one I measured in IDL
result = np.concatenate([result3,result2])


# In[ ]:


factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
empmass = factor*smaper/plxval**3
e_empmass = empmass*np.sqrt((esmaper/smaper)**2 +9.0*(plxprior/plxval)**2)


# In[ ]:


## this is mostly for checking things are reasonable
mka = kp - 5.0*(np.log10(1000.0/plxval)-1.)
mkb = ks - 5.0*(np.log10(1000.0/plxval)-1.)
a, b, c, d, e = result3
mka_err = 0.02## for now
mkb_err = 0.02## for now
mass1 = 10.0**(a + b*(mka-7.5) + c*(mka-7.5)**2 + d*(mka-7.5)**3 + e*(mka-7.5)**4)
mass2 = 10.0**(a + b*(mkb-7.5) + c*(mkb-7.5)**2 + d*(mkb-7.5)**3 + e*(mkb-7.5)**4)
mass1_err = (np.log(10)*(b+2*c*(mka-7.5)+3*d*(mka-7.5)**2+4*e*(mka-7.5)**3))*mass1*mka_err
mass2_err = (np.log(10)*(b+2*c*(mkb-7.5)+3*d*(mkb-7.5)**2+4*e*(mkb-7.5)**3))*mass2*mkb_err

model_err = np.sqrt(mass1_err**2+mass2_err**2)
model = mass1+mass2

mk1_tmp = mka[0] + mka_err*np.random.standard_normal(50000)#mka[0]+mka_err[0]*
mass1_tmp = a + b*(mk1_tmp-7.5) + c*(mk1_tmp-7.5)**2 + d*(mk1_tmp-7.5)**3 + e*(mk1_tmp-7.5)**4
mk2_tmp = mkb[0] + mkb_err*np.random.standard_normal(50000)#mka[0]+mka_err[0]*
mass2_tmp = a + b*(mk2_tmp-7.5) + c*(mk2_tmp-7.5)**2 + d*(mk2_tmp-7.5)**3 + e*(mk2_tmp-7.5)**4

#for i in range(0,len(empmass)):
#    print "{:10s}".format(name[i]), \
#    "{0:.3f}".format(empmass[i]),"{0:.3f}".format(e_empmass[i]), \
#    "{0:.4f}".format(model[i]),"{0:.4f}".format(model_err[i]),"{0:.3f}".format(100*model_err[i]/model[i]), \
#    "{0:.4f}".format(mka[i]),"{0:.4f}".format(mkb[i]), \
#    "{0:.3f}".format(ekp[i]),"{0:.3f}".format(eks[i]), \
#    "{0:.3f}".format(mass1[i]),"{0:.3f}".format(mass2[i]), \
#    "{0:.1f}".format(np.abs(empmass[i]-model[i])/np.sqrt(e_empmass[i]**2+model_err[i]**2))   


# In[ ]:


def lnlike(theta, smaper, esmaper, kp, ks, ekp, eks):
    au = 1.496e13
    msun = 1.989e33
    g = 6.6743e-8 
    a, b, c, d, e = theta[0:5]
    mplx = theta[5:theta.size]
    if np.min(mplx) <= 0:
        return -np.inf
    factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
    empmass = factor*smaper/mplx**3
    ## im not sure this is right with the plx prior in here
    ## seems like letting it float is the answer and not imposing the errors 
    ## a second time... try with 0.0 in front?
    ## wait, plxprior is not defined... is it 0?
    #e_empmass = empmass*np.sqrt((esmaper/smaper)**2 +9.*(plxprior/plxval)**2)
    e_empmass = empmass*(esmaper/smaper)**2
    mka = kp - 5.0*(np.log10(1000.0/mplx)-1.)
    mkb = ks - 5.0*(np.log10(1000.0/mplx)-1.)
    mass1 = 10.0**(a + b*mka + c*mka**2. + d*mka**3. + e*mka**4.)
    mass2 = 10.0**(a + b*mkb + c*mkb**2. + d*mkb**3. + e*mkb**4.)
    if np.min(mass1) <= 0:
        return -np.inf
    if np.min(mass2) <= 0:
        return -np.inf
    mka_err = ekp
    mkb_err = eks
    mass1_err = np.abs((np.log(10.)*(b+2.*c*mka+3.*d*mka**2+4.*e*mka**3.))*mass1*mka_err)
    mass2_err = np.abs((np.log(10.)*(b+2.*c*mkb+3.*d*mkb**2+4.*e*mkb**3.))*mass2*mkb_err)
    model_err = np.sqrt(mass1_err**2+mass2_err**2)
    model = mass1+mass2
    inv_sigma2 = 1.0/np.sqrt(e_empmass**2+model_err**2)
    return -0.5*(np.sum((empmass-model)**2*inv_sigma2 - np.log(inv_sigma2)))


# In[ ]:


def lnprior(theta, plxval, plxprior):
    mplx = theta[5:theta.size]
    lp = 0
    if np.min(mplx) <= 0:
        return -np.inf
    for i in range(0,len(mplx)):
        lp += ((np.float(mplx[i])-np.float(plxval[i]))**2)/(np.float(plxprior[i])**2)
    lp*=(-0.5)
    if not np.isfinite(lp):
        return 0.0
    return lp


# In[ ]:


def lnprob(theta, plxval, plxprior, smaper, esmaper, kp, ks, ekp, eks):
    lp = lnprior(theta, plxval, plxprior)
    if not np.isfinite(lp):
        return -np.inf
    like = lnlike(theta, smaper, esmaper, kp, ks, ekp, eks)
    if not np.isfinite(like):
        return -np.inf
    val = lp + like
    return val


# In[ ]:


ndim, nwalkers = result.size, 600
pos = [result + 5e-2*result*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, 
                                args=(plxval, plxprior, smaper, esmaper, kp, ks, ekp, eks),
                               threads=6)
## burn-in and/or testing
pos, prob, state = sampler.run_mcmc(pos, 5000)
sampler.reset()
print 'Finished burn/test phase'


# In[ ]:


import time
start_time = time.time()
nsteps = 500000
thin = 1000
kwargs = {'thin': thin }
print 'Starting run!'
for i, result in enumerate(sampler.sample(pos, iterations=nsteps, **kwargs)):
    if (i+1) % 5000 == 0:
        print("{0:5.1%}".format(float(i) / nsteps)),
        ("{0:5.2%}".format((time.time() - start_time)/60))
print 'Done, runtime:'
print (time.time() - start_time)/60
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))


# In[ ]:


import corner
from matplotlib.backends.backend_pdf import PdfPages
dat = sampler.chain
#print dat.shape,dat.size
#print dat.reshape(dat.size)
arr = dat.reshape(nwalkers*nsteps/thin,5+nbin)

fig = corner.corner(arr[:,0:5], labels=['a',r'b','c',r'd',r'e'], show_titles=True, title_kwargs={"fontsize": 11},title_fmt='.4f')
pp = PdfPages('output_params_log.pdf')
pp.savefig(fig)
pp.close()


# In[ ]:


dat = sampler.chain
prob = sampler.lnprobability
accept = sampler.acceptance_fraction
arr = dat.reshape((dat.shape)[0]*(dat.shape)[1],dat.shape[2])
for i in range(5,dat.shape[2]):
    print "{:10s}".format(name[i-5]),     "{0:.4f}".format(np.median(arr[:,i])),"{0:.4f}".format(np.std(arr[:,i])),"{0:.4f}".format(plxprior[i-5]),     "{0:.4f}".format(plxval[i-5]),"{0:.4f}".format((plxval[i-5]-np.median(arr[:,i]))/plxprior[i-5]),     "{0:.4f}".format((plxval[i-5]-np.median(arr[:,i]))/np.sqrt(plxprior[i-5]**2+np.std(arr[:,i])**2))
    


# In[ ]:


print a,b,c,d,e
factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
mass = factor*smaper/plxval**3
#e_empmass = empmass*np.sqrt((esmaper/smaper)**2 +9.0*(plxprior/plxval)**2)
#sma_au = sma*au*(1000/plxval)
#mass = (4.*(pi**2.))*((sma_au**3./per**2.)/g)/msun

mka = kp - 5.0*(np.log10(1000.0/plxval)-1.)
mkb = ks - 5.0*(np.log10(1000.0/plxval)-1.)

a = np.median(arr[:,0])
b = np.median(arr[:,1])
c = np.median(arr[:,2])
d = np.median(arr[:,3])
e = np.median(arr[:,4])
#mass1 = a + b*(mka-7.5) + c*(mka-7.5)**2 + d*(mka-7.5)**3 + e*(mka-7.5)**4
#mass2 = a + b*(mkb-7.5) + c*(mkb-7.5)**2 + d*(mkb-7.5)**3 + e*(mkb-7.5)**4
mass1 = 10.0**(a + b*mka + c*mka**2 + d*mka**3 + e*mka**4)
mass2 = 10.0**(a + b*mkb + c*mkb**2 + d*mkb**3 + e*mkb**4)
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


import matplotlib.pyplot as plt

rng = [np.min(np.concatenate([sysmass,sysmass_ben,mass])),
         np.max(np.concatenate([sysmass,sysmass_ben,mass]))]
plt.figure()
plt.plot(sysmass,mass,'ro')
plt.ylabel('Orbital Mass')
plt.xlabel('Predicted')
plt.plot(sysmass_ben,mass,'ro',color='b')
plt.plot(rng,rng)
plt.plot(sysmass_mann,mass,'ro',color='g')
plt.plot(rng,rng)
## blue = benedict
## green = mann idl
## red = python


# In[ ]:


#num = 7
#plt.plot(dat[2,:,num])
#print np.median(arr[:,num]),plxval[num-5]
newprob = prob.reshape((prob.shape)[0]*(dat.shape)[1])
print newprob.shape,arr[:,0].shape
n, bins, patches = plt.hist(newprob[np.isfinite(newprob)], 50, range=[10,45],normed=1, facecolor='green', alpha=0.75)
plt.show()


# In[ ]:


a = np.median(arr[np.isfinite(newprob),0])
b = np.median(arr[np.isfinite(newprob),1])
c = np.median(arr[np.isfinite(newprob),2])
d = np.median(arr[np.isfinite(newprob),3])
e = np.median(arr[np.isfinite(newprob),4])
print np.median(a),np.median(b),np.median(c),np.median(d),np.median(e)

rng = [np.min(np.concatenate([mka,mkb])),np.max(np.concatenate([mka,mkb]))]
print rng
mk = np.linspace(rng[0],rng[1],100)
mass1 = 10.0**(a + b*mka + c*mka**2 + d*mka**3 + e*mka**4)
a, b, c, d, e = result1
mass2 = a + b*(mk-7.5) + c*(mk-7.5)**2 + d*(mk-7.5)**3 + e*(mk-7.5)**4
mass3 = 0.585825902+3.87151019e-1*mk-1.21729937e-1*mk**2.+1.05529583e-2*mk**3.-2.72615872e-4*mk**4.

## red = new fit
## blue = Benedict
## green = How to constrain your M dwarf
plt.plot(mk,mass1,color='r')
plt.plot(mk,mass2,color='b')
plt.plot(mk,mass3,color='g')


for i in range(0,60):
    index = np.random.randint(len(arr[:,0]))
    if np.isfinite(newprob[index]):
        a, b, c, d, e = arr[index,0:5]
        #mass = a + b*(mk-7.5) + c*(mk-7.5)**2 + d*(mk-7.5)**3 + e*(mk-7.5)**4
        mass = 10.0**(a + b*mka + c*mka**2 + d*mka**3 + e*mka**4)
        plt.plot(mk,mass,color='r',lw=2,alpha=0.1)

plt.plot(mk,mass1,color='r')
plt.plot(mk,mass2,color='b')
plt.plot(mk,mass3,color='g')

plt.show()


# In[ ]:


## save out the relevant chains
import pyfits
pyfits.writeto('Mk-Mass_log_emcee.fits', sampler.chain, clobber=True)
pyfits.writeto('Mk-Mass_log_emcee_accept.fits', sampler.acceptance_fraction, clobber=True)
pyfits.writeto('Mk-Mass_log_emcee_lnprob.fits', sampler.lnprobability, clobber=True)
pyfits.writeto('Mk-Mass_log_emcee_acor.fits', sampler.acor, clobber=True)

