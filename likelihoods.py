import const
import background
import numpy as np
import SN
from scipy import interpolate

class DataBBN:

    def __init__(self,verbose=0):
        self.verbose = verbose
        self.Y_heobs = 0.2449
        self.Y_h2obs = 2.527
        self.sigma_heobs = 0.0040
        self.sigma_h2obs = 0.030
        self.sigma_h2th = 0.0060
        self.sigma_heth = 3.e-4

    def LnLike (self,BG):
        chi_yp = ( self.Y_heobs - BG.yp )**2 / ( self.sigma_heobs**2 + self.sigma_heth**2 )
        chi_H2 = ( self.Y_h2obs - (BG.x2p/1.e-5) )**2 / ( self.sigma_h2obs**2 + self.sigma_h2th**2 )
        if(self.verbose>0):
            print(' yp:%f' %BG.yp)
            print(' x2p:%f' %BG.x2p)
            print(' chi_yp:%f' %chi_yp)
            print(' chi_H2:%f' %chi_H2)
        return -0.5*(chi_yp+chi_H2)

class DataSNeIa:

    def __init__(self,verbose=0):
        self.verbose = verbose
        # Pantheon (alpha and beta not used - no nuisance parameters), fast
        self.like = SN.SN_likelihood(r'data/Pantheon/full_long.dataset')
        self.lnz1 = np.log(self.like.get_redshifts()+1)
        self.lnz1min = np.min(self.lnz1)
        self.lnz1max = np.max(self.lnz1)
        nbinnedlnz1 = 300
        self.binnedlnz1 = np.linspace(self.lnz1min,self.lnz1max,nbinnedlnz1)
        
    def LnLike(self,BG):
        
        def func(a): return BG.DeltaTau(a,1)*a*const.c/const.Mpc
        vfunc = np.vectorize(func,otypes=[float])
        
        binnedlnda = np.log(vfunc(np.exp(-self.binnedlnz1)))
        spl = interpolate.interp1d(self.binnedlnz1,binnedlnda,kind="cubic")
        return -self.like.loglike(np.exp(spl(self.lnz1)))

class DataBAO:

    def __init__(self,dataBAO,verbose=0):
        import inifile

        self.verbose = verbose

        Ini = inifile.IniFile("data/bao.dataset")
        self.dataBAO = dataBAO
        self.types = []
        self.zs = []
        self.Ds = []

        print('dataBAO:',self.dataBAO)

        i=0
        for section in self.dataBAO:
            self.types.append(Ini.ReadInt(section,"type"))
            self.zs.append(Ini.ReadFloat(section,"z"))
            if(self.types[-1]==1):
                self.Ds.append(Ini.ReadFloatArray(section,"DV"))
                i += 1
            elif(self.types[-1]==2):
                self.Ds.append(Ini.ReadFloatArray(section,"DM"))
                self.Ds.append(Ini.ReadFloatArray(section,"DH"))
                i += 2
        
#         # Alam et al. 2017 (arXiv:1607.03155)
#         self.nz0 = 3
#         self.nzobs0 = 2
#         self.rdrag0 = 147.78*const.Mpc
#         nobs = self.nz0*self.nzobs0
#         ncorr = nobs*(nobs+1)
#         self.z0 = np.array([0.38,0.51,0.61])
#         self.mean0 = np.array([1518,81.5,1977,90.4,2283,97.3])
#         err = np.array([22,1.9,27,1.9,32,2.1])
#         corr = np.array([1.0000,0.2280,1.0000,0.4970,0.1536,1.0000,0.1117,0.4873,0.2326,
#                          1.0000,0.1991,0.0984,0.5120,0.1571,1.0000,0.0520,0.2307,0.1211,0.5449,0.2408,1.0000])
#         covmat = np.empty([nobs,nobs])
#         for i in range(nobs):
#             for j in range(i+1):
#                 ij = i*(i+1)//2+j
#                 covmat[i,j] = err[i]*err[j]*corr[ij]
#                 covmat[j,i] = covmat[i,j]
#         self.invcov0 = np.linalg.inv(covmat)
#        
#         # Zarrouk et al. 2018 (arXiv:1801.03062)
#         self.nz1 = 1
#         self.nzobs1 = 2
#         nobs = self.nz1*self.nzobs1
#         ncorr = nobs*(nobs+1)
#         self.z1 = np.array([1.52])
#         self.mean1 = np.array([23.5e3,12.58])
#         err = np.array([1.78e3,0.68])
#         corr = np.array([1.0000,0.05,1.0000]) # sign of off diagonal term is flipped because Hr is inversely proportional to alpha_par
#         covmat = np.empty([nobs,nobs])
#         for i in range(nobs):
#             for j in range(i+1):
#                 ij = i*(i+1)//2+j
#                 covmat[i,j] = err[i]*err[j]*corr[ij]
#                 covmat[j,i] = covmat[i,j]
#         self.invcov1 = np.linalg.inv(covmat)
#        
#         # Bautista et al. 2017 (arXiv:1702.00176)
#         self.nz2 = 1
#         self.nzobs2 = 2
#         nobs = self.nz1*self.nzobs1
#         ncorr = nobs*(nobs+1)
#         self.z2 = np.array([2.33])
#         self.mean2 = np.array([9.07,37.77])
#         err = np.array([0.31,2.13])
#         corr = np.array([1.0000,0,1.0000]) # assumed diagonal
#         covmat = np.empty([nobs,nobs])
#         for i in range(nobs):
#             for j in range(i+1):
#                 ij = i*(i+1)//2+j
#                 covmat[i,j] = err[i]*err[j]*corr[ij]
#                 covmat[j,i] = covmat[i,j]
#         self.invcov2 = np.linalg.inv(covmat)

    def LnLike(self,BG):
        adrag = 1/(1+BG.zdrag)
        rdrag = BG.SoundHorizon(adrag)
        if(self.verbose>0):
            print(' z_drag:%f, r_drag[Mpc]:%f' %(BG.zdrag,rdrag/const.Mpc))
            
        lnl=0
        
        i=0
        for j in range(len(self.types)):
            type = self.types[j]
            z = self.zs[j]
            a = 1/(z+1)
            DM = BG.DeltaTau(a,1)*const.c
            DH = const.c*a*a*BG.dtauda(a)
            if(type==1):
                DV = (z*DM*DM*DH)**(1/3)
                x = (self.Ds[i][0]-DV/rdrag)/self.Ds[i][1]
                lnl += x*x
                i += 1
            elif(type==2):
                x = (self.Ds[i][0]-DM/rdrag)/self.Ds[i][1]
                lnl += x*x
                i +=1
                x = (self.Ds[i][0]-DH/rdrag)/self.Ds[i][1]
                lnl += x*x
                i += 1

#         f0 = rdrag/self.rdrag0
#         obs0 = np.empty(self.nz0*self.nzobs0)
#         for i in range(self.nz0):
#             obs0[2*i] = BG.DeltaTau(1/(1+self.z0[i]),1)/f0*const.c/const.Mpc # [Mpc]
#             obs0[2*i+1] = (1+self.z0[i])**2/BG.dtauda(1/(1+self.z0[i]))*f0*const.Mpc/const.km # in [km/sec/Mpc]
#         obs0 -= self.mean0
#         lnl += np.dot(obs0,np.dot(self.invcov0,obs0))
#
#         obs1 = np.empty(self.nz1*self.nzobs1)
#         for i in range(self.nz1):
#             obs1[2*i] = (1+self.z1[i])**2/BG.dtauda(1/(1+self.z1[i]))*rdrag/const.km # in [km/sec]
#             obs1[2*i+1] = BG.DeltaTau(1/(1+self.z1[i]),1)/rdrag*const.c/(1+self.z1[i])
#         obs1 -= self.mean1
#         lnl += np.dot(obs1,np.dot(self.invcov1,obs1))
#
#         obs2 = np.empty(self.nz2*self.nzobs2)
#         for i in range(self.nz2):
#             obs2[2*i] = (1+self.z2[i])**2/BG.dtauda(1/(1+self.z2[i]))*rdrag
#             obs2[2*i] = const.c/obs2[2*i]
#             obs2[2*i+1] = BG.DeltaTau(1/(1+self.z2[i]),1)/rdrag*const.c
#         obs2 -= self.mean2
#         lnl += np.dot(obs2,np.dot(self.invcov2,obs2))

        return -0.5*lnl

class DataH0:

    def __init__(self,verbose=0):
        self.verbose = verbose
        # Riess 2020 (arXiv:2001.03624)
        self.mean = 73.8
        self.error = 1.0

    def LnLike(self,BG):
        if(self.verbose>0):
            print(' H_0[km/sec/Mpc]:%f' %BG.H0())
        return -0.5*((BG.H0()-self.mean)/self.error)**2

class DataCMB:

    def __init__(self,verbose):
        self.verbose = verbose
        # model-independent prior on theta_MC according to Planck 2018 results VI (arXiv:1807.06209)
        self.mean = 0.010409
        self.error = 6e-6

    def LnLike(self,BG):
        astar = 1/(1+BG.zstar)
        rstar = BG.SoundHorizon(astar)
        rdiff = BG.SilkScale(astar)
        DA = BG.DeltaTau(astar,1)*const.c
        theta = rstar/DA
        if(self.verbose>0):
            print(' z_*:%f, r_*[Mpc]:%f'%(BG.zstar,rstar/const.Mpc))
            print(' 100*theta_*:%e'%(100*theta))
            print(' k_D[Mpc^{-1}]:%e'%(const.Mpc/rdiff))
            print(' 100*theta_diff:%e'%(100*rdiff*np.pi/DA))
        return -0.5*((theta-self.mean)/self.error)**2
    
class Likelihood:

    def __init__(self,useBAO,useH0,useCMB,useSNeIa,useBBN,dataBAO=None,verbose=0):
        self.useBAO = useBAO
        self.dataBAO = dataBAO
        self.useH0 = useH0
        self.useCMB = useCMB
        self.useSNeIa = useSNeIa
        self.useBBN = useBBN
        self.verbose = verbose

        self.nlikes = 1
        if(self.useBAO):
            self.BAO=DataBAO(self.dataBAO,self.verbose)
            self.nlikes += 1
        if(self.useH0):
            self.H0=DataH0(self.verbose)
            self.nlikes += 1
        if(self.useCMB):
            self.CMB=DataCMB(self.verbose)
            self.nlikes += 1
        if(self.useSNeIa):
            self.SNeIa=DataSNeIa(self.verbose)
            self.nlikes += 1
        if(self.useBBN):
            self.BBN=DataBBN(self.verbose)
            self.nlikes += 1

    def LnLike(self,BG,lnPprior):
        #import time

        #t0 = time.time()
        if(self.useBAO or self.useCMB):
            try: 
                BG.UpdateTherm()
            except ValueError:
                return [-np.inf for i in range(self.nlikes)]

        if(self.verbose>0):
            print("\n# likelihoods")
        lnL = [lnPprior]

        #t1 = time.time()
        if(self.useBAO):
            self.BAO.verbose = self.verbose
            lnL_BAO = self.BAO.LnLike(BG)
            lnL.append(lnL_BAO)
            lnL[0] += lnL_BAO

        #t2 = time.time()
        if(self.useH0):
            self.H0.verbose = self.verbose
            lnL_H0 = self.H0.LnLike(BG)
            lnL.append(lnL_H0)
            lnL[0] += lnL_H0

        #t3 = time.time()
        if(self.useCMB):
            self.CMB.verbose = self.verbose
            lnL_CMB = self.CMB.LnLike(BG)
            lnL.append(lnL_CMB)
            lnL[0] += lnL_CMB

        #t4 = time.time()
        if(self.useSNeIa):
            self.SNeIa.verbose = self.verbose
            lnL_SNeIa = self.SNeIa.LnLike(BG)
            lnL.append(lnL_SNeIa)
            lnL[0] += lnL_SNeIa

        #t5 = time.time()
        if(self.useBBN):
            self.BBN.verbose = self.verbose
            lnL_BBN = self.BBN.LnLike(BG)
            lnL.append(lnL_BBN)
            lnL[0] += lnL_BBN
            
        #t6 = time.time()
        #print(t5-t0,t1-t0,t2-t1,t3-t2,t4-t3,t5-t4,t6-t5)
        return lnL

