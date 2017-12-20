
# coding: utf-8

# In[1]:


import emcee
import matplotlib
import numpy as np
import corner
import scipy.optimize as op

#get_ipython().magic(u'matplotlib inline')


# In[2]:


# some important variables
g = 6.6743e-8
msun = 1.989e33
au = 1.496e13
pi = 3.14159


# In[3]:


#name = ['GJ1005','GJ22','GJ65','GJ1081','GJ234','G250-029','G193-027','GJ469','GJ473',
#       'GJ623','GJ748','GJ1245','GJ791.2','GJ831','Gl747','GJ2005']
#sma = np.array([0.3037,0.5106,2.0458,0.2712,1.0932,0.4417,0.1564,0.3139,0.9167,
#                0.2397,0.148,0.8267,0.1037,0.1448,0.2941,0.4619])
#esma = np.array([0.0005,0.0007,0.0066,0.0027,0.0007,0.0009,0.0008,0.0008,0.0017,
#                0.0014,0.0004,0.0008,0.0005,0.0005,0.0011,0.0006])
#per = np.array([1666.1,5694.2,9612,4066.1,6070.1,4946.3,1195.5,4223,5772.3,
#                  1367.4,901.7,6147,538.6,704.9,2102.151,17.3131*365.2422])
#eper = np.array([2.5,14.9,12,27.5,10.2,2.2,1.4,2.9,9.4,
#                 0.6,0.3,17,0.1,0.5,0.328,0.0097*365.2422])
#per *= 24*3600.0
#eper *= 24*3600.0
#k = np.array([6.39, 6.04, 5.33, 7.34, 5.49, 6.35, 7.78, 6.96, 6.04, 
#             5.92, 6.29, 6.85, 7.31, 6.38,  6.416,  8.241])
#delk = np.array([1.23, 1.93, 0.16, 0.96, 1.58, 1.03, 0.10, 1.01, 0.00, 
#                 2.87, 1.09, 1.05, 1.28, 1.20, 0.113,0.313])

## read in SMA, Period (Gaussian errors are fine?? confirm this)

## read in K magnitudes
#f = open('JHK_mags.txt','r')
#header1 = f.readline()
#i=0
#for line in f:
#    line = line.strip()
#    columns = line.split()
#    k[i] = float(columns[7])
#    ek[i] = float(columns[8])
#    i+=1
#f.close()
#ek = np.array([  0.02, 0.02, 0.02, 0.02, 0.04, 0.02, 0.02, 0.04, 0.02,
#                 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02])
#edelk = np.array([  0.02, 0.03, 0.02, 0.01, 0.02, 0.07, 0.02, 0.02, 0.01,
#                 0.14, 0.07, 0.01, 0.10, 0.01, 0.03, 0.01])

#plxval = np.array([166.6,99.2,373.7,65.2,241.0,95.6,112.6,76.4,235.5,
#                  125.0,98.4,219.9,113.4,125.3,122.3,128.49])
#plxprior = np.array([0.3,0.6,2.7,0.4,0.4,0.3,1.4,0.5,2.9,
#                    0.3,0.3,0.5,0.2,0.3,2.5,1.5])


# In[4]:


nbin = 50 ## number of binaries in sample

## new suggestion: put everything into ONE file within IDL, so we can just rin this K mag 
## script once and get everything in where we need it

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
result = np.concatenate([result1,result2])


# In[5]:


au = 1.496e13
msun = 1.989e33
g = 6.6743e-8 
factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
empmass = factor*smaper/plxval**3
e_empmass = empmass*np.sqrt((esmaper/smaper)**2 +9.0*(plxprior/plxval)**2)


# In[6]:


## this is mostly for checking things are reasonable
#empsma = sma*(1000.0/plxval)*au
#empmass = factor*smaper*1.0/plxval
#e_empmass = empmass*np.sqrt(3.0*(esma/empsma)**2. + 2.0*(eper/per)**2.)
#e_empmass = np.sqrt((esma*((12.0*np.pi**2*(empsma**2/per**2)/g)/msun))**2 + 
#            (eper*(6.0*np.pi**2*empsma**3/per**3)/g/msun)**2)
mka = kp - 5.0*(np.log10(1000.0/plxval)-1.)
mkb = ks - 5.0*(np.log10(1000.0/plxval)-1.)
a, b, c, d, e = result1
mass1 = a + b*(mka-7.5) + c*(mka-7.5)**2 + d*(mka-7.5)**3 + e*(mka-7.5)**4
mass2 = a + b*(mkb-7.5) + c*(mkb-7.5)**2 + d*(mkb-7.5)**3 + e*(mkb-7.5)**4
mka_err = 0.02## for now
mkb_err = 0.02## for now
mass1_err = np.sqrt((b*mka_err)**2+(2*c*mka_err*(mka-7.5))**2+(3*d*mka_err*(mka-7.5)**2)**2+
                    (4*e*mka_err*(mka-7.5)**3)**2)
mass2_err = np.sqrt((b*mkb_err)**2+(2*c*mkb_err*(mkb-7.5))**2+(3*d*mkb_err*(mkb-7.5)**2)**2+
                   (4*e*mkb_err*(mkb-7.5)**3)**2)    
model_err = np.sqrt(mass1_err**2+mass2_err**2)
model = mass1+mass2

mk1_tmp = mka[0] + mka_err*np.random.standard_normal(50000)#mka[0]+mka_err[0]*
mass1_tmp = a + b*(mk1_tmp-7.5) + c*(mk1_tmp-7.5)**2 + d*(mk1_tmp-7.5)**3 + e*(mk1_tmp-7.5)**4
mk2_tmp = mkb[0] + mkb_err*np.random.standard_normal(50000)#mka[0]+mka_err[0]*
mass2_tmp = a + b*(mk2_tmp-7.5) + c*(mk2_tmp-7.5)**2 + d*(mk2_tmp-7.5)**3 + e*(mk2_tmp-7.5)**4

for i in range(0,len(empmass)):
    print "{:10s}".format(name[i]),     "{0:.3f}".format(empmass[i]),"{0:.3f}".format(e_empmass[i]),     "{0:.4f}".format(model[i]),"{0:.4f}".format(model_err[i]),"{0:.3f}".format(100*model_err[i]/model[i]),     "{0:.4f}".format(mka[i]),"{0:.4f}".format(mkb[i]),     "{0:.3f}".format(ekp[i]),"{0:.3f}".format(eks[i]),     "{0:.3f}".format(mass1[i]),"{0:.3f}".format(mass2[i]),     "{0:.1f}".format(np.abs(empmass[i]-model[i])/np.sqrt(e_empmass[i]**2+model_err[i]**2))
    


# In[7]:


def lnlike(theta, smaper, esmaper, kp, ks, ekp, eks):
    au = 1.496e13
    msun = 1.989e33
    g = 6.6743e-8 
    a, b, c, d, e = theta[0:5]
    mplx = theta[5:theta.size]
    factor = (au**3.)*((4.0*np.pi**2.)/(g*msun))
    empmass = factor*smaper/plxval**3
    e_empmass = empmass*np.sqrt((esmaper/smaper)**2 +8.*(plxprior/plxval)**2)

    mka = kp - 5.0*(np.log10(1000.0/mplx)-1.)
    mkb = ks - 5.0*(np.log10(1000.0/mplx)-1.)
    mass1 = a + b*(mka-7.5) + c*(mka-7.5)**2 + d*(mka-7.5)**3 + e*(mka-7.5)**4
    mass2 = a + b*(mkb-7.5) + c*(mkb-7.5)**2 + d*(mkb-7.5)**3 + e*(mkb-7.5)**4
    mka_err = ekp
    mkb_err = eks
    mass1_err = np.sqrt((b*mka_err)**2+(2*c*mka_err*(mka-7.5))**2+(3*d*mka_err*(mka-7.5)**2)**2+
                       (4*e*mka_err*(mka-7.5)**3)**2)
    mass2_err = np.sqrt((b*mkb_err)**2+(2*c*mkb_err*(mkb-7.5))**2+(3*d*mkb_err*(mkb-7.5)**2)**2+
                       (4*e*mkb_err*(mkb-7.5)**3)**2)    
    model_err = np.sqrt(mass1_err**2+mass2_err**2)
    #print model_err[0],e_empmass[0],mka[0],mkb[0]
    if np.min(mass1) <= 0:
        return -np.inf
    if np.min(mass2) <= 0:
        return -np.inf
    model = mass1+mass2
    inv_sigma2 = 1.0/np.sqrt(e_empmass**2+model_err**2)
    #print empmass[0],model[0],e_empmass[0],inv_sigma2[0]
    return -0.5*(np.sum((empmass-model)**2*inv_sigma2 - np.log(inv_sigma2)))


# In[8]:


## test the lnlike formula (feed it Benedict and my own formula and compare likelihoods)
#plxval, plxprior, smaper, esmaper, kp, ks, ekp, eks
theta = result#[result1,plxval]
like1 = lnlike(theta, smaper, esmaper, kp, ks, ekp, eks)
tmp = np.array([0.2311,-0.1352, 0.0400, 0.0038, -0.0032]) # benedict fit value
tmpr = np.concatenate([tmp,result2])
theta = tmpr#[result1,plxval]
like2 = lnlike(theta, smaper, esmaper, kp, ks, ekp, eks)
print like1,like2


# In[9]:


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


# In[10]:


def lnprob(theta, plxval, plxprior, smaper, esmaper, kp, ks, ekp, eks):
    lp = lnprior(theta, plxval, plxprior)
    if not np.isfinite(lp):
        return -np.inf
    like = lnlike(theta, smaper, esmaper, kp, ks, ekp, eks)
    val = lp + like
    return val


# In[11]:


ndim, nwalkers = result.size, 200
pos = [result + 1e-2*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, 
                                args=(plxval, plxprior, smaper, esmaper, kp, ks, ekp, eks),
                               threads=3)
## burn-in and/or testing
pos, prob, state = sampler.run_mcmc(pos, 5000)
sampler.reset()


# In[ ]:


nsteps = 100000
for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
    if (i+1) % 5000 == 0:
        print("{0:5.1%}".format(float(i) / nsteps))
print 'Done!'


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
mass1 = a + b*(mka-7.5) + c*(mka-7.5)**2 + d*(mka-7.5)**3 + e*(mka-7.5)**4
mass2 = a + b*(mkb-7.5) + c*(mkb-7.5)**2 + d*(mkb-7.5)**3 + e*(mkb-7.5)**4
sysmass = mass1+mass2
a, b, c, d, e = result1
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
mass1 = a + b*(mk-7.5) + c*(mk-7.5)**2 + d*(mk-7.5)**3 + e*(mk-7.5)**4 ## new fit
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
        mass = a + b*(mk-7.5) + c*(mk-7.5)**2 + d*(mk-7.5)**3 + e*(mk-7.5)**4
        plt.plot(mk,mass,color='r',lw=2,alpha=0.1)
#    pl.plot(xl, m*xl+b, color="k", alpha=0.1)
#pl.plot(xl, m_true*xl+b_true, color="r", lw=2, alpha=0.8)
#pl.errorbar(x, y, yerr=yerr, fmt=".k")

plt.plot(mk,mass1,color='r')
plt.plot(mk,mass2,color='b')
plt.plot(mk,mass3,color='g')

plt.show()

