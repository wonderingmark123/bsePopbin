import warnings
import numpy as np
import matplotlib.pyplot as plt
import os 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.io as scio
import h5py
import math
import os
from tqdm import trange

def GetSFR(SFR,alpha,qav,binary,MminMINE,MminSFR):
    """
        Caculate the Binary Star forming rate of this program
        # Parameters

        SFR: float, original star forming rate obtained by obsercation.
        MminMINE: float or double, minimum higher initial mass.
        alpha: float or double, the exponent of mass function, which is different according to different version of mass function.
        qav: float, average mass ratio (q) of popbin program (default: 0.535)
        binary: int, Whether your SFR is for binary or all stars. ONLY 0 or 1!
        MminSFR: float, the minimum mass taken into account while caculating SFR.

        # Return

        Sb: float, Binary Star forming rate of this program
    """
    if(alpha<0):
        print('alpha should be larger than 0')
        Sb=None

    if(binary>0):
        ## whether single stars are counted in SFR or not 
        Sb=(MminMINE/MminSFR)**(1-alpha)*(alpha-2)*SFR /(1+qav) /MminSFR /(alpha-1)
    else:
        Sb=(MminMINE/MminSFR)**(1-alpha)*(alpha-2)*SFR /(3+qav) /MminSFR /(alpha-1)

    return(Sb)

# some basical parameters for caculation, the detailed information could be found in massFunction
alpha       = 2.7

## tmax=1d4*1d6; % in units of year

AllmassNUM  = 1e3
M1MAX       = 150
M1MIN       = 2
SFR         = GetSFR(43.5,alpha,0.535,0,M1MIN,5)
# SFR=10.6984
qNUM        = 93
aNUM        = 100
def massFunction(mass):
    """
        Caculate the birth rate of particular binary system.
        
        # Parameters
        ## Local

        mass: Double or float, the **higher** initial mass of binaries system.
        ## Global

        AllmassNUM: int, total mass number of Popbin.
        M1MAX: float or double, maximum higher initial mass 
        M1MIN: float or double, minimum higher initial mass 
        SFR: float or double, binary star forming rate 
        qNUM: int, total q (mass ratio) number of Popbin.
        aNUM: int, total a (seperation) number of Popbin.
        alpha: float or double, the exponent of mass function, which is different according to different version of mass function.

        # Return

        rate: float, the birth rate of particular binary system.
        
        # Warning 

        This function only support initial mass heavier than 1 solar mass! Otherwise the birth rate will be None!


    """ 
    global  AllmassNUM,M1MAX,M1MIN,SFR,qNUM,aNUM,alpha

    if (mass>1):
        rate=(mass/M1MIN)**(1-alpha)*SFR*(math.log(M1MAX)-math.log(M1MIN))/AllmassNUM/qNUM/aNUM*1e6 * (alpha-1);
        # rate=exp(-0.3.*mass);
    else:
        warnings.warn('Warning: This function only support initial mass heavier than 1 solar mass! Current mass is {}. The birth rate will be None!'.format(mass))
        rate=None

    return rate
def ImportOriData(
    NumParallel     = 1,
    DataFolder      = '..\Data',
    OutName         = '0415',
    SaveNPY         = False,
    LoadNPY         = False):
    """
    Import original data generated by Popbin_mine.f in ../Fortran0415, and delete some irrelated data. Finally output two array in the structured form.

    Please make sure your data files saved in ``{DataFolder}/{OutName}/all_Lx37erg_z0.01_rb_{OutName}.out``

    # Parameters
    NumParallel : int, the same as the setting in popbin program, (default: 1)
    DataFolder : string, the data saving folder of original output of popbin program.
    OutName : string, the same as the setting in popbin program, which is the version of popbin program.
    SaveNPY : bool, whether this program save the all_rb and all_wind array in NPY form in the same file folder. 
    LoadNPY : bool, whether this program load the all_rb and all_wind array in NPY form in the same file folder. Please set it to ``False`` in the first time. NPY form can speed up the loading progress significantly.

    # Returns

    RLselected: 2-D array, the meaning of each column are as folowing. The output array have been selected by ``DELallrbwind`` function.
    Windseleted: 2-D array, the meaning of each column are as folowing. The output array have been selected by ``DELallrbwind`` function.

    |Column| 0  |   1   |    2    |   3  |     4    |   5 |      6    |   7   |    8   |    9     |
    |  :----:  | :----:  | :----:  | :----:  | :----:  | :----:  | :----:  | :----:  | :----:  |:----:  | :----:  |
    |all_wind_1/WindSelected| m1|   m2|     Lx|     ecc0|   a|      tb0|    mt1|    mt2|    kw|   kw2|   
    |all_rb_1/RLselected| m1|   m2|     Lx|     ecc0|   a|      tb0|    mt1|    mt2|    kw|   kw2|    

    |Column| 10   |   11  |    12   |   13    |  14    |  15         |     16      |17    |  18     | 19|
    |  :----:  | :----:  | :----:  | :----:  | :----:  | :----:  | :----:  | :----:  | :----:  |:----:  | :----:  | 
    |all_wind_1/WindSelected| mx|     mx2|    eccx|   tbx|    f|      nnn|            uuu|    m_tr|   i   |    delta T|
    |all_rb_1/RLselected| mx|     mx2|    eccx|   tbx|    f|      m_transferrate| b_iso|  m_tr|   i   |  delta T|

    # Examples
    >>> from utils import *
    >>> all_rb_1,all_wind_1 = ImportOriData(
                            NumParallel     = 8,
                            DataFolder      = '..\Data',
                            OutName         = '0415',
                            SaveNPY         = True,
                            LoadNPY         = False)
    >>> all_rb_1,all_wind_1 = ImportOriData(
                            NumParallel     = 1,
                            DataFolder      = '..\Data',
                            OutName         = '0415',
                            SaveNPY         = False,
                            LoadNPY         = True)
    """
    if LoadNPY:
        RLselected = np.load(os.path.join(DataFolder,OutName,'RBdataSelected.npy'))
        Windseleted = np.load(os.path.join(DataFolder,OutName,'WinddataSelected.npy'))
        return RLselected[:,0:22],Windseleted
    RLori,Windori = np.array([]),np.array([])
    if NumParallel > 1:
        for i in trange(1,NumParallel+1):
            FileNameRL = os.path.join(DataFolder,OutName,str(i),
                ''.join([ 'all_Lx37erg_z0.01_rb_',OutName,'.out']))
            FilenameWind = os.path.join(DataFolder,OutName,str(i),
                ''.join([ 'all_Lx37erg_z0.01_wind_',OutName,'.out']))
            RLoriPart = np.loadtxt(FileNameRL)
            WindoriPart = np.loadtxt(FilenameWind)
            if np.size(RLori) == 0:
                RLori,Windori = RLoriPart,WindoriPart
            else:
                RLori = np.vstack((RLori,RLoriPart)).copy()
                Windori = np.vstack((Windori,WindoriPart)).copy()
        
    else:
        FileNameRL = os.path.join(DataFolder,
                ''.join([ 'all_Lx37erg_z0.01_rb_',OutName,'.out']))
        FilenameWind = os.path.join(DataFolder,
            ''.join([ 'all_Lx37erg_z0.01_wind_',OutName,'.out']))
        RLoriPart = np.loadtxt(FileNameRL)
        WindoriPart = np.loadtxt(FilenameWind)
        if np.size(RLori) == 0:
            RLori,Windori = RLoriPart,WindoriPart
    RLselected,Windseleted = DELallrbwind(RLori,Windori)
    if SaveNPY:
        np.save(os.path.join(DataFolder,OutName,'RBdataSelected'),RLselected)
        np.save(os.path.join(DataFolder,OutName,'WinddataSelected'),Windseleted)
    return RLselected[:,0:22],Windseleted
def DELallrbwind(all_rb, all_wind,
    minLx=39, 
    CHAbeamingFlag=False ,
    eddfac=1e4,
    DELtmax = 200,
    MaxDuration = 1):

    """
        Caculate the Binary Star forming rate of this program
        # Parameters

        all_rb: 2-D array, original data of Roche-lobe binaries obtained from Fortran0415 ``popbin_mine.f`` program.
        all_wind: 2-D array, original data of Wind Roche-lobe binaries obtained from Fortran0415 ``popbin_mine.f`` program.
        minLx: float or double, minimum luminosity of selection. (default: 39)
        CHAbeamingFlag=False ,
        eddfac=1e4,
        DELtmax = 200,
        MaxDuration = 1
        CHAbeamingFlag: bool, Whether changing beaming mode. ABANDONED! (default: False)
        eddfac: double, Maximum edding factor of selection (default: 1e4)
        DELtmax: float, Maximum evolution time of selection (default: 200 for Ring Galaxies)
        MaxDuration: float, Maximum stage duration of selection. In other words, the duration between two row of original bcm array in bse may be longer than several million year which is unrealiable for rapid mass transfer. (default: 1)


        # Return

        RLselected: 2-D array, Selected array
        Windseleted: 2-D array, Selected array

        # Criterion
        
        We delete following data which follows those criterion.
        - less than Maximum evolution time of selection
        - mass1 is smaller than 2Msun, which is unstable.
        - orbit is 0
        - $\delta T$ < 0
        - $\delta T$ > MaxDuration
        - Lx is lower than Lxmin

        # Warning 
        Please never set eddfac smaller than $10^4$ or set CHAbeamingFlag to True untill you understand this part of code!
    """
    # %% 1000Myr delete data whose max time is larger than DELtmax


    if (eddfac < 1e4):
        print(eddfac)

    if (CHAbeamingFlag):
        print(' beaming factor changed here')
    

    all_wind[:, 19] = all_wind[:, 7] - all_wind[:, 6]
    all_rb[:, 19] = all_rb[:, 7] - all_rb[:, 6]
    
    all_wind_1 = all_wind[all_wind[:, 7] <= DELtmax, :].copy()
    all_rb_1 =all_rb[all_rb[:, 7] <= DELtmax, :].copy()
    # %% change beaming factor
    if (CHAbeamingFlag):
        raise os.error('This part is not finished yet！')
        all_rb_1[:, 3] = all_rb_1[:, 3] * all_rb_1[:, 17]
        all_rb_1[:, 17] = getBeaming(2, all_rb_1[:, 16])
        all_rb_1[:, 3] = all_rb_1[:, 3] / all_rb_1[:, 17]
    

    # %% EDD factor
    # beaming_ori = all_rb_1(all_rb_1[:, 16] > eddfac, 17)
    # all_rb_1(all_rb_1[:, 16] > eddfac, 17) = getBeaming(1, eddfac)
    # all_rb_1(all_rb_1[:, 16] > eddfac, 3) = (all_rb_1(all_rb_1[:, 16] > eddfac, 3) 
    #     * beaming_ori / all_rb_1(all_rb_1[:, 16] > eddfac, 17))
    # % all_rb_1(all_rb_1(:,16)>eddfac,:)=[]
    # %% delete data whose mass1 is smaller than 2Msun , unstable

    all_rb_1 = all_rb_1[all_rb_1[:, 10] >= 2, :].copy()
    # % all_rb_1(all_rb_1(:,9)<2,:)=[]

    # %% sht try delete data whose  \delta T is small

    # % all_rb_1(all_rb_1(:,19)<1,:)=[]

    # %% orbit is 0
    # Condition = all_rb_1[:, 13]== 0
    all_rb_1 = all_rb_1[all_rb_1[:, 13]!= 0,:].copy()
    all_wind_1 = all_wind_1[all_wind_1[:, 13]!= 0, :].copy()
    # %% \delta T <0
    all_rb_1 = all_rb_1[all_rb_1[:, 19] > 0,:].copy()
    all_wind_1 = all_wind_1[all_wind_1[:, 19] > 0, :].copy()
    # delete \delta T > MaxDuration
    all_wind_1 = all_wind_1[all_wind_1[:, 19] <= MaxDuration, :].copy()
    all_rb_1 =all_rb_1[all_rb_1[:, 19] <= MaxDuration, :].copy()
    # %% Lx is lower than Lxmin
    all_rb_1 = all_rb_1[all_rb_1[:, 2] >= 10**minLx,:].copy()
    all_wind_1 = all_wind_1[all_wind_1[:, 2] >= 10**minLx, :].copy()
    # % all_rb_1(all_rb_1(:,9)<3,:)=[]
    # %% delete all mass ratio is larger than 3
    # % all_rb_1(all_rb_1(:,11)/all_rb_1(:,12)>3,:)=[]
    # % all_rb_1(all_rb_1(:,9)<2 & all_rb_1(:,16)>10,:)=[]
    # %% change the function of Lx
    # % all_rb_1(:,3)=0.1*(3.0d8)^2*all_rb_1(:,18)*Msun*1d7/yearsc
    # % all_wind_1(:,3)=0.1*(3.0d8)^2*all_wind_1(:,18)*Msun*1d7/yearsc

    # %% delete all eddfac larger than ???
    # % all_rb_1(all_rb_1(:,16)>=50,:)=[]
    return all_rb_1,all_wind_1



def Beamingfactor(beaming,b):
    """
        Different mode for beaming. The larger b is, the realistic the result is (maybe).

        # Parameter

        b: float, beaming factor.
        beaming: int, mode for caculation of beaming factor. The larger b is, the realistic the result is (maybe).(Defaut: 1)
        - 0 无修正：数目不做调整
        - 1 直接乘以b：数目乘以b=观测的数目
        - 2 方法二：像我昨天下午发的修正，乘以(8*b).^0.5/pi
        - 3 方法三，调整了一些估计项acos(1-b)*2/pi

        # Return 

        Beaming factor


    """
    if (beaming==0):
        xishu=1
    elif (beaming==1):
        xishu=b
    elif (beaming==2):
        # xishu = (8*b).^0.5/pi;
        xishu = math.sqrt( 2*b * ( 1 - b/2) )
    elif (beaming==3):
        xishu=math.acos(1-b)*2/math.pi;

    return xishu



def OBSplot(logFlag):
    """
        Plot the observation result of Wolter 2018 in seven Ring Galaxies.
        [2.6,2.5,2.3,20,8.0,4.1,4]
        |num|Ring Galaxies|SFR|
        |  :----:  | :----:  |
        |1|AM0644x3	|2.6|
        |2|Arp148x1	|2.5|
        |3|Arp143	|2.3|
        |4|Cartwheel	|20.0|
        |5|N922	|8.0|
        |6|Arp 147-Ring	|4.1|
        |7|Arp 284	|4.0|

        # Parameter

        logFlag: {bool,0,1}, whether return a log10(num) result or not. 
        obs: (local) 1-d array, the original observation data.

        # Return 

        This function return a turple
        obssort: array, the luminosity of each ULXs in seven Ring Galaxies in order.
        num: 1-D array, the number of ULXs more luminous than each luminosity of obssort.
        >>> num=np.linspace(len(obssort),1,len(obssort))
    """
    ## observation data
    obs=[[],[],[],[],[],[],[]];
    obs[0]=[8,12.6,28.6,9.8,22.4,2.1,2.3,5.1,8.5] #AM0644x3	
    obs[1]=[8.8,7.1,65.1] #Arp148x1	
    obs[2]=[2.7,9.8,2.9,1.9,3.8,1.6,4.3,1.4,1.2,2.8,7.5,1.4,2.2,1.6] #Arp143
    obs[3]=[3.36,11.86,11.15,65.49,7.43,6.37,17.7,2.66,15.93,19.47,2.3,1.19,2.47,1.45,0.86]
    ## Cartwheel
    obs[4]=[0.69,3.92,24.49,1.2,1.03,0.53,3.14,0.99,3.96,8.77,0.89,0.65]
    #N922
    obs[5]=[7.52,3.39,1.47,2.06,3.61,2.43,5.63,6.25,5.54] #Arp 147-Ring
    obs[6]=[64.1,0.5,0.6,0.7,0.8,0.6,0.3,0.8,24.1] #Arp 284
    obsall=[];
    SFRgalaxy= [2.6,2.5,2.3,20,8.0,4.1,4];

    ## analysis and plot
    for i in range (7):
        obsall=obsall+obs[i]
    obsall =np.log10( obsall) + 39
    obssort=np.sort(obsall);
    num=np.linspace(len(obssort),1,len(obssort))
    # num=log10(num);
    if(logFlag):
        num=np.log10(num);
    return (obssort,num)
    
def plot_Lx_mesaNS(beaming,
    tmax = 200,
    Number_XLFpoint = 50,
    minLx = 39,
    maxLx = 48,
    fileFolder = '../Data/PaperPlotData/'):
    """
        Different mode for beaming. The larger b is, the realistic the result is (maybe).

        # Parameter

        beaming: float, beaming factor.
        tmax: float, maximum evolution time for bianries.
        Number_XLFpoint: int, the length of the array while ploting XLF.
        minLx: float, minimum luminosity of XLF
        maxLx: float, maximum luminosity of XLF
        fileFolder: path, the saving path of mesa data 

        # Return 

        NNsumwind: array, the number of predicted ULXs via Roche-Lobe overflow by mesa program which is more luminous than particular luminosity.
        NNsumrb: array, the number of predicted ULXs via Wind Roche-Lobe overflow by mesa program which is more luminous than particular luminosity.

        # Warning

        make sure that there is Lx_N1207.npy and hanchen_first_neutron1124.out in your file folder!
    """
    ## beaming_lower_limit=1d-3;
    ApplyNSEddington = False
    Lx_jian=(maxLx-minLx)/Number_XLFpoint;
    Lx_x=np.linspace(minLx,maxLx,Number_XLFpoint)
    lx_sun=3.845*1e33

    all = np.loadtxt(os.path.join(fileFolder,'hanchen_first_neutron1124.out')) # 11556行19列

    if(ApplyNSEddington):
        mesa_lx = np.loadtxt('D:/study/bsegit/bse-pap/Lx_N_m_Porb.log')
    else:
        mesa_lx = np.load(os.path.join(fileFolder,'Lx_N1207.npy')) # 900774行10列
        # mesa_lx = np.loadtxt('./Lx_N1207.log') # 900774行10列
        mesa_lx[:, 3:6]= mesa_lx[:, 3:6] *lx_sun

    ## fake all_rb_1
    all_rb_1=np.zeros((len(mesa_lx),20))
    all_rb_1[:,2]   = mesa_lx[:,3];
    all_rb_1[:,16]  = mesa_lx[:,6];
    all_rb_1[:,13]  = mesa_lx[:,8];
    all_rb_1[:,10]  = mesa_lx[:,9];
    all_rb_1[:,6]   = mesa_lx[:,7];
    all_rb_1[:,19]  = mesa_lx[:,2];
    # plot_Lx_b_mesa(all_rb_1,100,100);

    ## delete ultraluminous data
    ## delete evolution time largerr than tmax
    delnum=[]
    for i in range (len(mesa_lx)):
        if (mesa_lx[i,6]<1e-2) or (mesa_lx[i,7] / 1e6 > tmax):
            delnum.append(i)
    mesa_lx=np.delete(mesa_lx,delnum,axis=0)      
    del(delnum); del(i);

    ## note
    # (massI)( obJ ( delt  (Lapp_total)   Lapp_overflow Lapp_wind)   beamingfactor
    # 1       2       3        4             5                 6       7
    # all=load('D:\study\bsegit\MW\hmxb.dat');

    # 1  2   3  4   5  6   7   8   9  10 11 12  13    14 15 16  17  18 19
    # m1,m2,Lx,ecc0,a,tb0,mt1,mt2,kw,kw2,mx,mx2,eccx,tbx,f,nnn,uuu,m_tr,i
    # m1,m2,Lx,ecc0,a,tb0,mt1,mt2,kw,kw2,mx,mx2,eccx,tbx,f,m_transferrate,b_iso,m_tr,i
    # 20 \delta T
    aNUM=100;
    qNUM=93;


    # SFR=10.6984;

    yearsc=3.1557e07;
    Msun=1.9891e30;
    yearday=365.0;

    all_1=all;
    m_loc=4;
    tb_loc=6;

    m_loc=11;
    tb_loc=14;

    # mline=[9.625:mstep:20.625,21.25:0.5:30.25];
    mline = np.arange(2-0.125,29.5+0.125+0.01,0.25)
    tbline=np.arange(0.2 - 0.05,3.0 + 0.05+0.01,0.1)
    # mline=[0.1:mstep:20.625];
    # tbline=(0.001:0.1:3.0);

    # change orbit  to circle
    ecc = all_1[:, 12]
    all_1[:,tb_loc-1] = all_1[:,tb_loc-1] *yearday * (1 - ecc**2 )** 1.5

    delnum=[]
    for i in range (len(all_1)):
        if (all_1[i,m_loc-1]<=min(mline)) or (all_1[i,m_loc-1]>=max(mline)) or (all_1[i,tb_loc-1]<= 10**(min(tbline))) or (all_1[i,tb_loc-1]>= 10**(max(tbline))):
            delnum.append(i)

    all_1=np.delete(all_1,delnum,axis=0)
    del(delnum); del(i);

    mpos=[]
    for i in range(len(all_1)):
        for j in range(len(mline)):
            if all_1[i,m_loc-1]<mline[j]:
                mpos.append(j)
                break
    tbpos=[]
    for i in range(len(all_1)):
        for j in range(len(tbline)):
            if all_1[i,tb_loc-1]<10**tbline[j]:
                tbpos.append(j)
                break

    mlen=len(all_1);
    # result=zeros(length(mline),length(tbline));
    RateMorb=np.zeros([len(tbline),len(mline)]);
    for i in range (mlen):
        RateMorb[tbpos[i]-1,mpos[i]-1]=RateMorb[tbpos[i]-1,mpos[i]-1]+massFunction(all_1[i,0])


    RateMorb = RateMorb[0:np.shape(RateMorb)[0]-1,:][:,0:np.shape(RateMorb)[1]-1]
    ## wind 
    NN=np.zeros(Number_XLFpoint);
    data=mesa_lx;

    mlen=len(mesa_lx);
    for i in range (mlen):
        wei=math.ceil((math.log10(data[i,5]) - minLx) / Lx_jian);
        if(wei<1):
            continue
        
        ##     deal with wei which is out of range coz the wrong chosen number of
        #     maxLx
        if(wei>len(NN)):
            wei=len(NN)

        massI= data[i,0];
        obJ= data[i,1];
    #   NN(wei)=NN(wei)+data(i,20)*massFunction(data(i,1)).*Beamingfactor(beaming,data(i,17));
    
        NN[wei-1]=NN[wei-1]+data[i,2]*RateMorb[int(obJ)-1, int(massI)-1]

    NNwind = NN;


    ## rl 
    NN=np.zeros(Number_XLFpoint);
    data=mesa_lx;
    mlen=len(mesa_lx);
    for i in range (mlen):
        wei=math.ceil((math.log10(data[i,4]) - minLx) / Lx_jian);
        if (wei<1):
            continue
    ##     deal with wei which is out of range coz the wrong chosen number of
    #     maxLx
        if(wei>len(NN)):

            wei=len(NN)

        massI= data[i,0]
        obJ= data[i,1]
        #   NN(wei)=NN(wei)+data(i,20)*massFunction(data(i,1)).*Beamingfactor(beaming,data(i,17));
    
        NN[wei-1]=NN[wei-1]+data[i,2]*RateMorb[int(obJ)-1, int(massI)-1]*Beamingfactor(beaming,data[i,6])
    NNrb = NN;
    NNsumrb = (np.cumsum(NNrb[::-1]))[::-1];

    NNsumwind = (np.cumsum(NNwind[::-1]))[::-1];

    ## check 
    # hist(log10(data(:, 7)))

    ## lx-beaming
    return NNsumwind,NNsumrb

def plot_Lx_b_mesa(data,xge,yge):
    # ABANDONED
    ml = 3;
    m2 = 11;
    beaming=2;
    # 0 无修正：数目不做调整
    # 1 直接乘以b：数目乘以b=观测的数目
    # 2 方法二：像我昨天下午发的修正，乘以(8*b).^0.5/pi
    # 3 方法三，调整了一些估计项acos(1-b)*2/pi

    data[:,ml-1]=math.log10(data[:,ml-1]);
    massline=np.linspace(min(data[:,ml-1]),max(data[:,ml-1]),xge);
    orbline=np.linspace(min(math.log10(data[:,m2-1])),max(math.log10(data[:,m2-1])),yge);
    NN=np.meshgrid(orbline,massline)*0;
    mlen=len(data);
    for i in range (mlen):
    
        masswei=np.ceil((data[i,ml-1]-massline[1])/(massline[2]-massline[1]));
        orbwei=np.ceil((math.log10(data[i,m2-1])-orbline[1])/(orbline[2]-orbline[1]));
        if(masswei==0):
            masswei=1
        if(orbwei ==0):
            orbwei=1
        if(data[i,14]>100):
            #NN(masswei,orbwei)=NN(masswei,orbwei)+data(i,20)...
            #*massFunction(data(i,1)).*Beamingfactor(beaming,data(i,17));
            NN[masswei-1,orbwei-1]=NN[masswei-1,orbwei-1]+data[i,19]*Beamingfactor(beaming,data[i,16]);
        else:
            #NN(masswei,orbwei)=NN(masswei,orbwei)+...
            #data(i,20)*massFunction(data(i,1));
            NN[masswei-1,orbwei-1]=NN[masswei-1,orbwei-1]+data[i,19];

    NN[NN<np.max(NN)/1000]=0;

    # a=surf(massline,orbline,log10(NN'));
      # a.FaceColor = 'interp';
    # color_mine=fliplr(gray')';
    # colormap(color_mine)
    # a.EdgeColor='none';
    # set(a,'AlphaData',NN==0)
    # view(2)
    # ylabel('log10(Luminosity)/log10(erg/s)')
    # xlabel('log10(beaming)')
    # colorbar('LineWidth',1.5);
      # a=pcolor(orbline,massline,log10(NN));
      # a=pcolor(massline,orbline,log10(NN'));
      # xlim([0.3 1.3])
      # ylim([-1 4])
    # box on

    ax = plt.subplot(111, projection='3d')
    ax.set_title('???');
    ax.plot_surface(massline,orbline,math.log10(NN),rstride=2, cstride=1, cmap=plt.cm.Spectral)
    #设置坐标轴标签
    ax.set_xlabel('log10(beaming)')
    ax.set_ylabel('log10(Luminosity)/log10(erg/s)')
    #ax.set_zlabel('C')
    plt.show()

    return NN
    
def plot_Lx_NNsum_gedian(data,beaming,
    Number_XLFpoint = 50,
    minLx = 39,
    maxLx = 48):
    """
        Generate data for ploting the X-ray Luminosity function(XLF). 

        # Parameter

        data: 2-D array, in the same form with all_rb_1 or all_wind_1, and could be part of them.
        beaming: float, beaming factor.
        Number_XLFpoint: int, the length of the array while ploting XLF.
        minLx: float, minimum luminosity of XLF
        maxLx: float, maximum luminosity of XLF

        # Return 
        (NNsum,Lx_x,np.log10(NNsum))
        NNsum: 1-D array, the number of predicted ULXs which is more luminous than particular luminosity (Lx_x) obtained from the input data.
        Lx_x: 1-D array, particular luminosity (Lx_x) 
        >>> Lx_x=np.linspace(minLx,maxLx,Number_XLFpoint)

    """
    Lx_jian=(maxLx-minLx)/Number_XLFpoint;
    Lx_x=np.linspace(minLx,maxLx,Number_XLFpoint);
    NN=np.zeros(Number_XLFpoint);
    mlen=len(data);
    for i in range (mlen):
        wei=math.ceil((math.log10(data[i,2])-minLx)/Lx_jian);
        if(wei<0):
            continue
        ##  deal with wei which is out of range coz the wrong chosen number of  maxLx
        if(wei>len(NN)-1):
            wei=len(NN);
            #    continue;
        #     NN(wei)=NN(wei)+data(i,20)*massFunction(data(i,1)).*Beamingfactor(beaming,data(i,17));
        if(data[i,14]>100):
            NN[wei-1]=NN[wei-1]+data[i,19]*massFunction(data[i,0])*Beamingfactor(beaming,data[i,16])
        else:
            NN[wei-1]=NN[wei-1]+data[i,19]*massFunction(data[i,0])

    NNsum=(np.cumsum(NN[::-1]))[::-1]

    # plot(Lx_x,log10(NN))

    # after cumsum of NN. (more smooth [\doge] )

    #plt.figure()
    #plt.plot(Lx_x,np.log10(NNsum))
    #plt.xlabel('log10(Lx)/(erg/s)')
    #plt.ylabel('log10(N)')
    #plt.show()

    return (NNsum,Lx_x,np.log10(NNsum))

def plot_lx_tb_gedian_paper(data_ori, Number_Lx, Number_orb, my_title, Limit_cri=1000):
    """
        Generate data for ploting the distribution of luminosity-orbit plane.

        # Parameter

        data_ori: 2-D array, in the same form with all_rb_1 or all_wind_1, and could be part of them.
        beaming: float, beaming factor.
        Number_Lx: int, the length of the array related to luminosity.
        Number_orb: int, the length of the array related to orbit.
        my_title: string, title of this figure
        Limit_cri: number, the minimum distribution for showing in the plot. For example, if Limit_cri=1000 and the maximum distribution is 1. then the plot will not display distribution less than 1/1000.

        # Return 
        NN: 2-D array, the distribution number of each panel.

    """

    beaming=1;
    # 0 无修正：数目不做调整
    # 1 直接乘以b：数目乘以b=观测的数目
    # 2 方法二：像我昨天下午发的修正，乘以(8*b).^0.5/pi
    # 3 方法三，调整了一些估计项acos(1-b)*2/pi
    ml=3
    data = data_ori.copy()
    data[:,ml-1]=np.log10(data[:,ml-1])
    Lxline=np.linspace(38.99,41,Number_Lx)
    orbline=np.linspace(-1.62,4.225,Number_orb)
    NN=np.zeros((len(Lxline),len(orbline)))
    mlen=len(data)
    for i in range (mlen):
        masswei=np.ceil((data[i,ml-1]-Lxline[0])/(Lxline[1]-Lxline[0]))
        orbwei=np.ceil((np.log10(data[i,13])-orbline[0])/(orbline[1]-orbline[0]))
        if(masswei > Number_Lx or orbwei > Number_orb or masswei<=0 or orbwei <= 0):
            continue
        if(data[i,14]>100):
            NN[int(masswei)-1,int(orbwei)-1]=NN[int(masswei)-1,int(orbwei)-1]+data[i,19]*massFunction(data[i,0])*Beamingfactor(beaming,data[i,16])
        else:
            NN[int(masswei)-1,int(orbwei)-1]=NN[int(masswei)-1,int(orbwei)-1]+data[i,19]*massFunction(data[i,0])
    
    cri=(np.max(NN))/Limit_cri
    for i in range (np.shape(NN)[0]):
        for j in range (np.shape(NN)[1]):
            if NN[i,j]<cri:
                NN[i,j]=0

    Lx_m,orb_m=np.meshgrid(Lxline,orbline)
    NN_m=np.log10(np.transpose(NN))

    plt.figure()
    vmax=np.ceil(np.max(NN_m[np.isfinite(NN_m)])*2)/2
    vmin=np.floor(np.min(NN_m[np.isfinite(NN_m)])*2)/2
    norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
    h=plt.pcolor(Lx_m,orb_m,NN_m,cmap=plt.cm.YlOrRd,norm=norm)
    plt.colorbar(h,ticks=np.linspace(vmin,vmax,int((vmax-vmin)*2+1)))
    plt.xlabel('log10(Luminosity) / log10(erg/s)')
    plt.ylabel('log10($P_{orb}$)/log10(day)')
    plt.title(my_title)

    return NN
def plot_mass_tb_gedian_paper(data_ori, Number_mass, Number_orb, my_title, Limit_cri=1000):
    """
        Generate data for ploting the distribution of mass-orbit plane.

        # Parameter

        data_ori: 2-D array, in the same form with all_rb_1 or all_wind_1, and could be part of them.
        beaming: float, beaming factor.
        Number_mass: int, the length of the array related to mass.
        Number_orb: int, the length of the array related to orbit.
        my_title: string, title of this figure
        Limit_cri: number, the minimum distribution for showing in the plot. For example, if Limit_cri=1000 and the maximum distribution is 1. then the plot will not display distribution less than 1/1000.

        # Return 
        NN: 2-D array, the distribution number of each panel.

    """
    ml=11;

    beaming=1;
    # 0 无修正：数目不做调整
    # 1 直接乘以b：数目乘以b=观测的数目
    # 2 方法二：像我昨天下午发的修正，乘以(8*b).^0.5/pi
    # 3 方法三，调整了一些估计项acos(1-b)*2/pi
    data = data_ori.copy()
    # data[:,ml-1]=np.log10(data[:,ml-1]);
    massline=np.linspace(min(data[:,ml-1])-np.spacing(1.0)*100,max(data[:,ml-1]),Number_mass)
    orbline=np.linspace(min(np.log10(data[:,13]))-np.spacing(1.0)*100,max(np.log10(data[:,13])),Number_orb)
    NN=np.zeros((len(massline),len(orbline)))
    mlen=len(data)
    for i in range(mlen):
        masswei=np.ceil((data[i,ml-1]-massline[0])/(massline[1]-massline[0]))
        orbwei=np.ceil((np.log10(data[i,13])-orbline[0])/(orbline[1]-orbline[0]))
        if (data[i,14]>100):
            xishu=Beamingfactor(beaming,data[i,16])
            rate=massFunction(data[i,0])*xishu
            NN[int(masswei)-1,int(orbwei)-1] = NN[int(masswei)-1,int(orbwei)-1]+data[i,19]*rate
        else:
            NN[int(masswei)-1,int(orbwei)-1]=NN[int(masswei)-1,int(orbwei)-1]+data[i,19]*massFunction(data[i,0])
    
    cri=(np.max(NN))/Limit_cri
    for i in range (np.shape(NN)[0]):
        for j in range (np.shape(NN)[1]):
            if NN[i,j]<cri:
                NN[i,j]=0

    mass_m,orb_m=np.meshgrid(massline,orbline)
    NN_m=np.log10(np.transpose(NN))
    plt.figure()

    vmax=np.ceil(np.max(NN_m[np.isfinite(NN_m)])*2)/2
    vmin=np.floor(np.min(NN_m[np.isfinite(NN_m)])*2)/2
    norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
    h=plt.pcolor(mass_m,orb_m,NN_m,cmap=plt.cm.YlOrRd,norm=norm)
    plt.colorbar(h,ticks=np.linspace(vmin,vmax,int((vmax-vmin)*2+1)))
    plt.xlabel('log10($M_{Donor}$) / $M_{solar}$')
    plt.ylabel('log10($P_{orb}$)/log10(day)')
    plt.title(my_title)
    #plt.xlim([])
    #plt.ylim([0,3])
    #plt.title('time of duration(MW)')
    #plt.savefig('4-23M-P(MW)1')


    #a=surf(massline,orbline,log10(NN'));
    #color_mine=fliplr(gray')';
    #colormap(color_mine)
    #a.EdgeColor='none';
    #set(a,'AlphaData',NN==0)
    #view(2)
    #ylabel('log10(P_{orb})/log10(day)')
    #xlabel('log10(M_{Donor}) / M_{solar}')
    #colorbar('LineWidth',1.5);
    #box on

    return NN

def SETplot(picName):
    # ABANDONED

    # title(picName);
    #set(gca,'FontName','Times New Roman','FontSize',25,'LineWidth',1.5);
    #set(gca, 'XMinorTick', 'on')
    #set(gca, 'YMinorTick', 'on')
    #print(gcf, '-djpeg', '-opengl', sprintf('-r%d',600), [picName,'.jpg']);
    pass
