import sys
import numpy as np

ddir = sys.argv[1]
name = sys.argv[2]#'GJ1210AB'
keck = sys.argv[3]#'_Nep=13_dT=10.06'
nsteps = int(sys.argv[4])#100000
nwalkers = int(sys.argv[5])#200
thin = int(sys.argv[6])#100

#print type(nsteps),type(nwalkers),type(thin)

#name = 'GL417B'
#keck = '_Nep=12_dT=13.23'

#name = 'GL569B'
#keck = '_Nep=19_dT=11.87'

dtor = np.pi/180.0

#
# read in binary sep and PA data
#
dat1 = np.genfromtxt(ddir+name+keck+'_sep-PA.txt',usecols=(0,1,2,3,4),skip_header=1)
ep1 = dat1[:,[0]]
se  = dat1[:,[1]]
pa  = dat1[:,[2]] * dtor
ese = dat1[:,[3]]
epa = dat1[:,[4]] * dtor

#
# determine initial guess
#
guess = np.genfromtxt(ddir+name+keck+'_guess.txt',usecols=(0),skip_header=1)
# convert angles from degrees to radians
guess[[3,4,5]] = guess[[3,4,5]] * dtor
# convert from t_0 to lambda_ref
tref = 2455197.5
guess[6] = 2.0 * np.pi * (tref - guess[6]) / guess[0] + guess[5]
# ensure angles are between 0 and 2pi
guess[[4,5,6]] = (guess[[4,5,6]] + 2.0*np.pi) % (2.0*np.pi)
asc_guess = guess[4]
arg_guess = guess[5]
guess[4] = arg_guess - asc_guess
guess[5] = arg_guess + asc_guess 

#
# define likelihood function for visual binary parameters
#
def lnlike(theta, se, pa, ese, epa, ep1, tref):
    per, sma, ecc, inc, omo, opo, lam = theta
    arg = (opo + omo)*0.5
    asc = (opo - omo)*0.5
    # calculate mean anomaly
    MA1 = 2.0 * np.pi * (ep1 - tref) / per + lam - arg
    # iteratively solve for eccentric anomaly
    niter = 4
    EA1 = MA1 + ecc * np.sin(MA1) + 0.5 * ecc**2 * np.sin(2.0*MA1)
    for i in range(1,niter):
        EA1 += (MA1 - EA1 + ecc * np.sin(EA1))/(1.0 - ecc * np.cos(EA1))
    # calculate the Thiele-Innes orbital elements
    sinarg = np.sin(arg)
    cosarg = np.cos(arg)
    sinasc = np.sin(asc)
    cosasc = np.cos(asc)
    cosinc = np.cos(inc)
    A =  cosarg*cosasc - sinarg*sinasc*cosinc
    B =  cosarg*sinasc + sinarg*cosasc*cosinc
    F = -sinarg*cosasc - cosarg*sinasc*cosinc
    G = -sinarg*sinasc + cosarg*cosasc*cosinc
    X1 = np.cos(EA1) - ecc
    Y1 = np.sin(EA1) * np.sqrt(1.0 - ecc**2)
    # calculate Cartesian orbit positions (dimensionless)
    dr1 = B*X1 + G*Y1
    dd1 = A*X1 + F*Y1
    # compute orbit separation (dimensionless) and position angle
    se1 = sma*np.sqrt(dr1**2 + dd1**2)
    pa1 = np.arctan2(dr1,dd1)
    # compute log likelihood
    dse = se-se1
    dpa = np.arctan2(np.sin(pa-pa1),np.cos(pa-pa1))
    inv_sigma2 = 1.0/np.concatenate([ese,epa])**2
    return -0.5*(np.sum(np.concatenate([dse,dpa])**2*inv_sigma2 - np.log(inv_sigma2)))

#
# define priors on orbital parameters
#
def lnprior(theta):
    per, sma, ecc, inc, omo, opo, lam = theta
    arg = (opo + omo)*0.5
    asc = (opo - omo)*0.5
    if per <= 0.0:
        return -np.inf
    if sma <= 0.0:
        return -np.inf
    if ecc < 0.0 or ecc >= 0.999:
        return -np.inf
    if inc <= 0.0 or inc >= np.pi:
        return -np.inf
    if asc < -np.pi or asc >= 3.0*np.pi:
        return -np.inf
    if arg < -np.pi or arg >= 3.0*np.pi:
        return -np.inf
    if lam < -np.pi or lam >= 3.0*np.pi:
        return -np.inf
    return np.log(1.0/per * 1.0/sma * np.sin(inc))

#
# define final log likelihood function
#
def lnprob(theta, se, pa, ese, epa, ep1, tref):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, se, pa, ese, epa, ep1, tref)


#
# initialize and run sampler
#
ndim = 7
nwalkers = nwalkers
scat = [ 1e-4*guess[0], 1e-4*guess[1], 1e-4, 1e-4, 1e-4, 1e-4, 1e-4 ]
pos = [guess + scat*np.random.randn(ndim) for i in range(nwalkers)]

print "prior on guess (check that this is not -inf)", lnprior(guess)
import emcee
kwargs = {'thin': thin }
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(se, pa, ese, epa, ep1, tref), threads=8)
import time
start_time = time.time()
sampler.run_mcmc(pos, nsteps, **kwargs)
print time.time() - start_time
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

    #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(data, nplanets, priorstruc, inposindex, struc1, args), threads=nthreads)
    #sampler.run_mcmc(pos, nsteps)

#
# save chain
#
import pyfits
pyfits.writeto(ddir+name+keck+'_emcee.fits', sampler.chain, clobber=True)
pyfits.writeto(ddir+name+keck+'_emcee_accept.fits', sampler.acceptance_fraction, clobber=True)
pyfits.writeto(ddir+name+keck+'_emcee_lnprob.fits', sampler.lnprobability, clobber=True)
pyfits.writeto(ddir+name+keck+'_emcee_acor.fits', sampler.acor, clobber=True)


#print sampler.chain.shape

#
# make a fancy corner plot
#
import corner
from matplotlib.backends.backend_pdf import PdfPages

dat = sampler.chain
radeg = 57.295776
arr = dat.reshape(nwalkers*nsteps/thin,7)
arr[:,0]=arr[:,0]/365.2422
arr[:,3]=arr[:,3]*radeg
omo=arr[:,4]*radeg
opo=arr[:,5]*radeg
lam=arr[:,6]*radeg
arg = (opo + omo)*0.5
asc = (opo - omo)*0.5
arr[:,4] = arg
arr[:,5] = asc

fig = corner.corner(arr[:,0:5], labels=['Period (yr)',r'$a$ (arcsec)','e',r'inc (deg)',r'arg (deg)',r'asc (deg)'], show_titles=True, title_kwargs={"fontsize": 11},title_fmt='.4f')
pp = PdfPages(name+keck+'.pdf')
pp.savefig(fig)
pp.close()
