import sys
import const
import background
import numpy as np
import likelihoods
import inifile
import mcmc

def main():
    args = sys.argv
    if(len(args)<2):
        print("error: the number of input parameters is not correct; input command must be")
        print(">$ python inifiles.py input_file_name")
        sys.exit(1)

    Ini = inifile.IniFile(args[1])
    Ini.Dump()

    section = "OUTPUT"
    root = Ini.ReadString(section,"root")
    Ini.Dump(root+"_params.ini")
    
    # model calculation with fiducial parameters
    section = "COSMOLOGY"
    BG = background.Background(Ini.ReadInt(section,"neutrino_hierarchy"),Ini.ReadInt(section,"wtype"),verbose=1)
    paramsfid = np.array([Ini.ReadFloat(section,"obh2"),Ini.ReadFloat(section,"odmh2"),Ini.ReadFloat(section,"odeh2"),
                          Ini.ReadFloat(section,"nnu"),Ini.ReadFloat(section,"mnu"),Ini.ReadFloat(section,"w[0]"),
                          Ini.ReadFloat(section,"w[1]"),Ini.ReadFloat(section,"w[2]"),Ini.ReadFloat(section,"w[3]")])
    BG.SetParams(paramsfid)

    # test
    BG.OutputComovingDistance(root+'_distance.txt',1e-3,1e3,100)
    
    # likelihood calculation
    section = "LIKELIHOODS"
    LF = likelihoods.Likelihood(Ini.ReadBoolean(section,"use_BAO"),Ini.ReadBoolean(section,"use_H0"),
                                Ini.ReadBoolean(section,"use_CMB"),Ini.ReadBoolean(section,"use_SNeIa"),
                                Ini.ReadBoolean(section,"use_BBN"),
                                dataBAO=Ini.ReadString(section,"data_BAO").split(),verbose=1)
    lnL = LF.LnLike(BG,0)
    print(" ln(L)=",lnL[:])

    # MCMC
    inp = input("\nDoes everything seem going fine? Then let's run MCMC [Y/n]:")
    if(inp=='n'):
        sys.exit()
    section ='MCMC'
    BG.verbose  = LF.verbose  = 0
    paramrange = [Ini.ReadFloatArray(section,"obh2"),Ini.ReadFloatArray(section,"odmh2"),Ini.ReadFloatArray(section,"odeh2"),
                  Ini.ReadFloatArray(section,"nnu"),Ini.ReadFloatArray(section,"mnu"),Ini.ReadFloatArray(section,"w[0]"),
                  Ini.ReadFloatArray(section,"w[1]"),Ini.ReadFloatArray(section,"w[2]"),Ini.ReadFloatArray(section,"w[3]")]
    map_varied = np.array([i for i in range(len(paramrange)) if len(paramrange[i]>0)])
    range_varied = np.array([paramrange[i] for i in map_varied])
    driver = mcmc.MCMC(paramsfid,Ini.ReadBoolean(section,"parallel"))
    driver.SetParams(map_varied,range_varied,Ini.ReadInt(section,"nwalkers"),verbose=1)
    chain = root+".h5"
    driver.Run(chain,Ini.ReadInt(section,"nsteps"),BG,LF)
    
if __name__ == '__main__':
    main()
