import sys
import background
import numpy as np
import likelihoods
import emcee
import time

def LnPost(params_varied,*args):
    if(np.all(params_varied[:]>=args[0].rangev[:,0]) and np.all(params_varied[:]<=args[0].rangev[:,1])):
        p = np.array(args[0].pf)
        p[args[0].mapv[:]] = params_varied[:]
        args[1].SetParams(p)
        lnp = np.sum(args[0].rangev[:,2]*np.log(np.abs(params_varied[:]))) # prior
        blobs = args[2].LnLike(args[1],lnp)
        blobs.extend(args[1].GetDerivedParams())
    else:
        blobs = [-np.inf for i in range(args[0].nblobs+1-args[1].ndparams)]
        blobs.extend([ np.nan for i in range(args[1].ndparams)])

    return (*blobs,)

class MCMC():
    import time
    
    def __init__(self,params_fid,parallel=True,nprocesses=1):
        from multiprocessing import cpu_count
        self.pf = np.array(params_fid)
        self.parallel = parallel
        if(self.parallel):
            self.nprocesses = min(cpu_count(),nprocesses) if(nprocesses>0) else cpu_count()
    
    def SetParams(self,map_varied,range_varied,nwalkers,verbose=0):
        self.nv = len(map_varied)
        self.mapv = np.array(map_varied)
        self.rangev = np.array(range_varied[:,:])
        self.nwalkers = nwalkers
        self.verbose = verbose
        if(self.verbose>0):
            print("\n# MCMC")
            print(" number of walkers:",self.nwalkers)
        
    def Run(self,chainfilename,nsteps,BG,LF):
        backend = emcee.backends.HDFBackend(chainfilename)
        try: # restart
            print(" number of existing steps:",backend.iteration)
            print(" restarting from existing chains")
        except: # new chain
            print(" creating brand-new chains")
            backend.reset(self.nwalkers,self.nv)
        
        # initial condition
        p0 = np.random.randn(self.nwalkers,self.nv)
        p0[:,:] = p0[:,:]*self.rangev[None,:,3] +self.pf[None,self.mapv[:]]
        
        # metadata (derived parameters etc.)
        dtype = []
        if(LF.useBAO):
            dtype.append(("lnL_BAO",float))
        if(LF.useH0):
            dtype.append(("lnL_H0",float))
        if(LF.useCMB):
            dtype.append(("lnL_CMB",float))
        if(LF.useSNeIa):
            dtype.append(("lnL_SNeIa",float))
        if(LF.useBBN):
            dtype.append(("lnL_BBN",float))
        BG.SetDerivedParams(LF.useBAO or LF.useCMB)
        for i in range(BG.ndparams):
            dtype.append((BG.dparams[i],float))
        self.nblobs = len(dtype)

        if(self.parallel):
            from multiprocessing import Pool
            import os
            
            # following classes are required to be created so that LnPost should be pickable
            BG0 = background.Background(BG.nu.hierarchy,BG.de.wtype,verbose=0)
            BG0.SetDerivedParams(BG.with_therm)
            LF0 = likelihoods.Likelihood(LF.useBAO,LF.useH0,LF.useCMB,LF.useSNeIa,LF.useBBN,LF.dataBAO,verbose=0)
            MC0 = MCMC(self.pf)
            MC0.SetParams(self.mapv,self.rangev,self.nwalkers,verbose=0)
            MC0.nblobs = self.nblobs
                
            # Intel MKL should be serialized as instructed in emcee userguide
            os.environ["OMP_NUM_THREADS"] = "1"
            os.environ["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"
            
            with Pool(processes=self.nprocesses) as pool:
                
                self.sampler = emcee.EnsembleSampler(MC0.nwalkers,MC0.nv,LnPost,args=(MC0,BG0,LF0),backend=backend,blobs_dtype=dtype,pool=pool)
                print(" parallel run")
                print(" {0} processes".format(self.nprocesses))
                self.sampler.run_mcmc(p0,nsteps,progress=True)
        else:
            
            print(" serial run")
            self.sampler = emcee.EnsembleSampler(self.nwalkers,self.nv,LnPost,args=(self,BG,LF),blobs_dtype=dtype)
            self.sampler.run_mcmc(p0,nsteps,progress=True)

