import sys
import numpy as np
from scipy import integrate,interpolate,optimize
import const

class BBN:
    
    def __init__(self):
        filename = "data/bbn.dat" # based on PArthENoPE taken from CLASS 
        data = np.loadtxt(filename)
        self.spl_yp = interpolate.interp2d(data[:,0],data[:,1],data[:,2],kind="linear") # cubic doesn't work; why?
        
    def yp(self,obh2,dnnu):
        return self.spl_yp(obh2,dnnu)

class MassiveNu:

    def __init__(self,hierarchy):
        self.nmass = 3 # number of mass eigen states
        self.nnu = const.nnu_standard
        self.mass = np.zeros(self.nmass)
        self.hierarchy = hierarchy

        self.lammax = 1e3 # max of am/T
        self.lammin = 1e-3 # min of am/T
        
        nlam = 1000
        lnlam = np.linspace(np.log(self.lammin),np.log(self.lammax),nlam)
        lnrho = np.zeros(nlam)
        for i in range(nlam):
            lam = np.exp(lnlam[i])
            lnrho[i] = integrate.quad(lambda x:x**2*np.sqrt(x**2+lam**2)/(np.exp(x)+1),0,100)[0]
        lnrho = np.log(lnrho[:]*const.nu_energy) # ratio of massive to massless
        self.spl_lnrho = interpolate.make_interp_spline(lnlam,lnrho)
        
    def Rho(self,a):
        lam = a*self.mass[:]*const.eV/const.TCnuB
        r = np.empty(self.nmass)
        for i in range(self.nmass):
            if(lam[i]>self.lammax): # non-relativistic
                r[i] = lam[i]*const.nu_number
            elif(lam[i]<self.lammin):
                r[i] = 1
            else:
                r[i] = np.exp(self.spl_lnrho(np.log(lam[i])))
        return r
        
    def SetParams(self,nnu,mnu): # mnu [eV]
        self.nnu = nnu
        mlower = 0
        mupper = mnu
        mlight = optimize.brentq(lambda x:self.lightest2tot(x)-mnu,mlower,mupper)
        if(self.hierarchy==1): # normal
            self.mass[0] = mlight
            self.mass[1] = np.sqrt(mlight**2+const.m2nu21)
            self.mass[2] = np.sqrt(mlight**2+const.m2nu21+const.m2nu32)
        elif(self.hierarchy==-1): # inverted
            self.mass[0] = np.sqrt(mlight**2-const.m2nu21+const.m2nu32)
            self.mass[1] = np.sqrt(mlight**2+const.m2nu32)
            self.mass[2] = mlight
        else: # degenerate
            self.mass[0:2] = mlight
        
    def lightest2tot(self,mlight):
        if(self.hierarchy==1): # normal
            return mlight+np.sqrt(mlight**2+const.m2nu21)+np.sqrt(mlight**2+const.m2nu21+const.m2nu32)
        elif(self.hierarchy==-1): # inverted
            return np.sqrt(mlight**2-const.m2nu21+const.m2nu32)+np.sqrt(mlight**2+const.m2nu32)+mlight
        else: # degenerate
            return mlight*self.nmass

class DarkEnergy:
    
    def __init__(self,wtype):
        self.wtype = wtype

        self.lnrinterp = False
        self.binned = False
        if(self.wtype==0):
            pass
        else:
            self.lnrinterp = True

            if(self.wtype>=10): # binned w
                self.binned = True
                self.pbin = np.floor(self.wtype/10)
                self.nbin = self.wtype%10+1
            else:
                print("error: not supported yet")
                sys.exit(1)
                
            self.lnz1min = 0
            self.lnz1max = np.log(3e3)
            if(self.binned):
                self.nlnz1 = self.nbin+1
            else:
                self.nlnz1 = 50

        self.w0 = -1
        self.wa = 0
        
    def SetParams(self,wparams):
        self.wparams = np.array(wparams)

        if(self.lnrinterp):
            if(self.binned):
                self.abin = np.linspace(np.exp(-self.lnz1min),np.exp(-self.lnz1max),self.nlnz1)**(1/self.pbin)
                lnz1 = -np.log(self.abin)
            else:
                lnz1 = np.linspace(self.lnz1min,self.lnz1max,self.nlnz1)
                
            lnr = np.empty(self.nlnz1)
            lnr[0] = 0
            for i in range(1,self.nlnz1):
                lnr[i] = lnr[i-1]+3*integrate.quad(lambda x:1+self.EoS(np.exp(-x)),lnz1[i-1],lnz1[i])[0]
            self.spl = interpolate.make_interp_spline(lnz1,lnr)
        
        # effective w0 and wa
        if(self.lnrinterp):
            def func(lnz1,w0,wa) :return 3*((1+w0+wa)*lnz1-wa*(1-np.exp(-lnz1)))
            opt,cov = optimize.curve_fit(func,lnz1,lnr)
            self.w0 = opt[0]
            self.wa = opt[1]
        elif(self.wtype==0):
            self.w0 = self.wparams[0]
            self.wa = self.wparams[1]
            
    def DumpParams(self):
        print(" wtype:",self.wtype)
        print(" binned w?:",self.binned)
        if(self.binned):
            print(" binned z:",1/self.abin-1)
        print(" wparams:",self.wparams)
        print(" interpolation?:",self.lnrinterp)
            
    def EoS(self,a):
        if(self.wtype==0):
            return self.wparams[0]+(1-a)*self.wparams[1]
        elif(self.binned):
            i = int((1-a**self.pbin)*self.nbin)
            return self.wparams[i]
        
    def Rho(self,a):
        if(self.lnrinterp):
            return np.exp(self.spl(-np.log(a)))
        elif(self.wtype==0):
            return np.exp(-3*((1+self.wparams[0]+self.wparams[1])*np.log(a)+self.wparams[1]*(1-a)))
            
class Background:
    
    def __init__(self,hierarchy,wtype,verbose=0):
        self.verbose = verbose
        # default cosmological parameters
        self.ogh2 = np.pi**2/15*const.TCMB**4/(const.c*const.hbar)**3/const.rhoch2
        self.nu = MassiveNu(hierarchy)
        self.onuh2_nomass = self.ogh2*7/8*(const.TCnuB/const.TCMB)**4*self.nu.nnu
        self.obh2 = 0.0224
        self.odmh2 = 0.120
        self.odeh2 = 0.311
        self.de = DarkEnergy(wtype)
        self.bbn = BBN()
        self.yp = 0.24714 # based on PRIMAT assuming obh2 = 0.02236, neff = 3.046
        
        # derived parameters
        self.SetDerivedParams()
        
    def SetParams(self,params):
        #params is an array containing [obh2,odmh2,odeh2,Gamma_Gyr,mratio,nnu,mnu]
        self.obh2 = params[0]
        self.odmh2 = params[1]
        self.odeh2 = params[2]
        self.nu.SetParams(params[3],params[4])
        self.onuh2_nomass = self.ogh2*7/8*(const.TCnuB/const.TCMB)**4*self.nu.nnu
        self.de.SetParams(params[5:])
        self.yp = self.bbn.yp(self.obh2,self.nu.nnu-const.nnu_standard)[0]
        
        if(self.verbose>0):
            print("\n# cosmological parameters:")
            print(" omega_g*h^2: %e"%self.ogh2)
            print(" omega_nu*h^2 (no mass): %e"%self.onuh2_nomass)
            print(" omega_dm*h^2 (no decay): %e"%self.odmh2)
            print(" omega_de*h^2: %e"%self.odeh2)
            print(" neutrino mass hierarchy [1:NH,-1:IH],0:NH]:",self.nu.hierarchy)
            print(" neutrino masses [eV]:",self.nu.mass[:])
            print(" neutrino NR redshifts:",const.TCnuB/const.eV/self.nu.mass[:])
            print(" yp: %e"%self.yp)
            self.de.DumpParams()
            
    def dtauda(self,a):
        x = 1/np.sqrt(self.ogh2+self.onuh2_nomass*np.average(self.nu.Rho(a))+(self.obh2+self.odmh2)*a+self.odeh2*a**4*self.de.Rho(a))
        return x/const.BigH
        
    def H0(self):
        return 100/self.dtauda(1)/const.BigH
        
    def DeltaTau(self,a1,a2):
        return integrate.quad(self.dtauda,a1,a2)[0]

    def drsda(self,a):
        cs = np.sqrt(1/3/(1+0.75*a*self.obh2/self.ogh2))*const.c
        return cs*self.dtauda(a)
    
    def SoundHorizon(self,a):
        return integrate.quad(self.drsda,0,a)[0]
    
    def UpdateTherm(self):
        from HyRec import pyrec
        
        if(self.verbose>0):
            print("thermal history is computed by HyRec")
        
        ok = 0
        pyrec.rec_build_history_wrap(const.TCMB/const.kB,self.obh2,self.odmh2,ok,self.odeh2,self.de.w0,self.de.wa,self.yp,self.nu.nnu,tuple(self.nu.mass))
        
        # initial guess based on Hu & Sugiyama 1996
        obh2 = self.obh2
        om = self.obh2+self.odmh2
        g1 = 0.0783*obh2**-0.238/(1+39.5*obh2**0.763)
        g2 = 0.560/(1+21.1*obh2**1.81)
        self.zstar = 1048*(1+0.00124*obh2**-0.738)*(1+g1*om**g2)
        b1 = 0.313*om**-0.419*(1+0.607*om**0.674)
        b2 = 0.238*om**0.223
        self.zdrag = 1345*om**0.251/(1+0.659*om**0.828)*(1+b1*obh2**b2)
        
        zstar1 = self.zstar*0.95
        zstar2 = self.zstar*1.05
        zdrag1 = self.zdrag*0.95
        zdrag2 = self.zdrag*1.05
        
        akthom = obh2*const.rhoch2/const.c*(1-self.yp)/const.m_H*const.sigmaT
        
        nthrm = 10
        zlower = min(zstar1,zdrag1)
        zupper = max(zstar2,zdrag2)
        zthrm = np.linspace(zlower,zupper,nthrm)
        tau_opt = np.empty(nthrm)
        tau_drag = np.empty(nthrm)
        tau_opt[0] = integrate.quad(lambda z:
                    pyrec.hyrec_xe(1/(z+1))*akthom*self.dtauda(1/(z+1)),0,zthrm[0])[0]
        tau_drag[0] = +integrate.quad(lambda z:
                    pyrec.hyrec_xe(1/(z+1))*akthom*self.dtauda(1/(z+1))/(0.75*self.obh2/self.ogh2/(1+z)),0,zthrm[0])[0]
        for i in range(1,nthrm):
            tau_opt[i] = tau_opt[i-1]+integrate.quad(lambda z:
                                pyrec.hyrec_xe(1/(z+1))*akthom*self.dtauda(1/(z+1)),zthrm[i-1],zthrm[i])[0]
            tau_drag[i] = tau_drag[i-1]+integrate.quad(lambda z:
                                pyrec.hyrec_xe(1/(z+1))*akthom*self.dtauda(1/(z+1))/(0.75*self.obh2/self.ogh2/(1+z)),zthrm[i-1],zthrm[i])[0]
        spl_opt = interpolate.make_interp_spline(np.log(zthrm[1:]),np.log(tau_opt[1:]))
        spl_drag = interpolate.make_interp_spline(np.log(zthrm[1:]),np.log(tau_drag[1:]))
        
        self.zstar = optimize.brentq(lambda x:spl_opt(np.log(x)),zstar1,zstar2)
        self.zdrag = optimize.brentq(lambda x:spl_drag(np.log(x)),zdrag1,zdrag2)
        
    def SetDerivedParams(self):
        self.ndparams = 2
        self.dparams = ["H0","Age"]
        
    def GetDerivedParams(self):
        H0 = self.H0()
        t_in_Gyr = integrate.quad(lambda x:x*self.dtauda(x),1e-3,1)[0]/const.Gyr
        return [H0,t_in_Gyr]
