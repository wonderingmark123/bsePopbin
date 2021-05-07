import warnings
from utils import *
import numpy as np
import os
from paper_plot import PlotAll

def main0415():
    print('Loading data now')
    all_rb_1,all_wind_1 = ImportOriData(
    NumParallel     = 8,
    DataFolder      = '..\Data',
    OutName         = '0415',
    SaveNPY         = True,
    LoadNPY         = False)
    PlotAll(all_rb_1,all_wind_1,
    FileFolder = 'D:/study/bsegit/bsePopbin/Data/PaperPlotData')
    
    return

def main1117():
    PlotAll([],[],
    FileFolder = 'D:/study/bsegit/bsePopbin/Data/PaperPlotData')
    return
def Examples():
    print('Loading data now')
    all_rb_1,all_wind_1 = ImportOriData(
    NumParallel     = 8,
    DataFolder      = '..\Data',
    OutName         = '0415',
    SaveNPY         = True,
    LoadNPY         = False)

    allinall=np.vstack((all_rb_1, all_wind_1))


    # Plot mesa luminosity function 
    [N_NSwind,N_NSrb] = plot_Lx_mesaNS(1,fileFolder='D:/study/bsegit/bsePopbin/Data/PaperPlotData');
    Lx_x=np.linspace(39,48,50);
    plt.figure()
    plt.plot( Lx_x , np.log10(N_NSwind),label='NS mesa')

    # Plot observation luminosity function of star-burst galaxies
    plt.figure()
    (obssort,num)=OBSplot(1)
    plt.scatter(obssort,num,marker='s',color=(0.6350,0.0780,0.1840),label='Obsevation')


    # plot XLF of BH WRLOF
    plt.figure()
    temparray=all_wind_1[all_wind_1[:,9]==14,:].copy()
    N_BHwind    = plot_Lx_NNsum_gedian(temparray,0)[0]

    # Plot mass-orb plane of BH RL
    plt.figure()
    temparray=all_rb_1[all_rb_1[:,9]==14,:].copy()
    plot_mass_tb_gedian_paper(temparray,200,100,' M2-Porb-BH-RLOF');

    # Plot lx- tb plane of NS WRLOF
    plt.figure()
    temparray=all_wind_1[all_wind_1[:,9]==13,:].copy()
    plot_lx_tb_gedian_paper(temparray,200,100,' Lx-Porb-NS-wind')
    return
if __name__ == '__main__':
    main0415()