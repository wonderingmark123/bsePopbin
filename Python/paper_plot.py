import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.io as scio
import h5py
import math
import os

from utils import *
def PlotAll(all_rb_1,all_wind_1,
    FileFolder = 'D:/study/bsegit/bsePopbin/Data/PaperPlotData'):
    # 字体预设1
    font1 = {'family':'Times New Roman','weight':'normal','size':16,}
    if np.size(all_rb_1) == 0:

        

        ###########读取数据############
        mat = h5py.File(os.path.join(FileFolder,'1117_paper.mat'))
        # mat文件里可能有多个cell，各对应着一个dataset
        # 可以用keys方法查看cell的名字
        # print(mat.keys())

        # 0     1       2       3       4       5       6       7       8       9         
        # m1,   m2,     Lx,     ecc0,   a,      tb0,    mt1,    mt2,    kw,   kw2,    
        # m1,   m2,     Lx,     ecc0,   a,      tb0,    mt1,    mt2,    kw,   kw2,    
        # 10      11      12      13      14      15              16      17      18      19
        # mx,     mx2,    eccx,   tbx,    f,      nnn,            uuu,    m_tr,   i       delta T
        # mx,     mx2,    eccx,   tbx,    f,      m_transferrate, b_iso,  m_tr,   i       delta T

        all_rb_1=mat['all_rb_1']
        all_rb_1=np.transpose(all_rb_1)
        all_wind_1=mat['all_wind_1']
        all_wind_1=np.transpose(all_wind_1)

    ########## main ##########

    ##### initial parameters

    yearsc      =3.1557e07
    Msun=1.9891e30
    minLx=39
    maxLx=48
    # xi=Mdot/(sum(sum(data(:,1:2)))/tmax)/tmax*1d6
    CHAbeamingFlag=0
    eddfac=1e4
    changeBeta = 1

    ##### all in all. ONE for ALL!!!

    all_rb_1 = all_rb_1[all_rb_1[:,10]>0,:];
    all_wind_1 = all_wind_1[all_wind_1[:,10]>0,:];

    allinall=np.vstack((all_rb_1, all_wind_1))

    allinall = allinall[allinall[:,10]>0,:];
    ##### Luminosity function

    ge=50
    beaming=1

    Lx_jian=(maxLx-minLx)/ge;
    Lx_x=np.linspace(minLx,maxLx,ge);

    ## different maximum evolution time for binaries

    tspan = [5,10,15,20,50,100,200];

    N_time=np.zeros((len(tspan),50))
    for i in range(len(tspan)):
        tmax = tspan[i];
        [N_NSwind, N_NSrb] = plot_Lx_mesaNS(beaming,tmax,fileFolder=FileFolder)

        
        temparray=all_rb_1[(all_rb_1[:,6]<=tmax) + (all_rb_1[:,9]==14),:].copy()
        N_BHrb = plot_Lx_NNsum_gedian(temparray,beaming)[0]


        temparray=all_wind_1[(all_wind_1[:,6]<=tmax) + (all_wind_1[:,9]==14),:].copy()
        N_BHwind    = plot_Lx_NNsum_gedian(temparray,0)[0]

        N_time[i]=N_BHrb + N_BHwind + N_NSrb + N_NSwind

    plt.figure()
    (obssort,num)=OBSplot(1)
    plt.scatter(obssort,num,marker='s',color=(0.6350,0.0780,0.1840),label='Obsevation')
    plt.plot(Lx_x, np.log10(N_time[2]),linestyle='-',color=(0.3010,0.7450,0.9330),label='15Myr')
    plt.plot(Lx_x, np.log10(N_time[3]),linestyle='-.',color=(0.4660,0.6740,0.1880),label='20Myr')
    plt.plot(Lx_x, np.log10(N_time[4]),linestyle=':',color=(0.3010,0.250,0.9),label='50Myr')
    plt.plot(Lx_x, np.log10(N_time[5]),linestyle='--',color=(0.9,0.2,0.6),label='100Myr')
    plt.plot(Lx_x, np.log10(N_time[6]),linestyle='-',color=(0,0.4470,0.7410),label='200Myr')
    plt.xlabel('log10(Luminosity) / log10(erg/s)',font1)
    plt.ylabel('log10(Number of ULX)',font1)
    plt.legend(prop=font1)
    plt.yscale('linear')
    plt.xlim((39,40.5))
    plt.ylim((0.5,2))

    #plt.show() #记得注释掉！！！！
    ##

    [N_NSwind,N_NSrb] = plot_Lx_mesaNS(beaming,fileFolder=FileFolder);

    temparray=all_rb_1[all_rb_1[i,9]==14,:].copy()
    N_BHrb      = plot_Lx_NNsum_gedian(temparray,beaming)[0]

    selectnum=[]
    for i in range (len(all_wind_1)):
        if all_wind_1[i,9]==14:
            selectnum.append(i)
    temparray=all_wind_1[selectnum,:]   
    del(selectnum); del(i);
    N_BHwind    = plot_Lx_NNsum_gedian(temparray,0)[0]

    N_all       = N_BHrb + N_BHwind + N_NSrb + N_NSwind

    plt.figure()
    plt.plot( Lx_x , np.log10(N_all),linestyle='-',color=(0,0.4470,0.7410),label='Total BPS')
    (obssort,num)=OBSplot(1)
    plt.scatter(obssort,num,marker='s',color=(0.6350,0.0780,0.1840),label='Obsevation')
    plt.plot(Lx_x,np.log10( N_BHrb),linestyle='--',color=(0.9010,0.550,0.2),label='RL BH')
    plt.plot(Lx_x, np.log10(N_NSrb),linestyle='--',color=(0.4660,0.6740,0.1880),label='RL NS')
    plt.plot(Lx_x, np.log10(N_BHwind),linestyle=':',color=(0.9010,0.550,0.2),label='Wind BH')
    plt.xlabel('log10(Luminosity) / log10(erg/s)',font1)
    plt.ylabel('log10(Number of ULX)',font1)
    plt.legend(prop=font1)
    plt.yscale('linear')
    plt.xlim((39,40.5))
    plt.ylim((0.5,2))
    #plt.LegendBox='on';
    #plt.FileName="MAINLuminosityFunction.jpg";
    #plt.show() #########
    a=1


    plt.figure()
    plt.plot( Lx_x , np.log10(N_all),linestyle='-',color=(0,0.4470,0.7410),label='Total BPS')
    (obssort,num)=OBSplot(1);
    plt.scatter(obssort,num,marker='s',color=(0.6350,0.0780,0.1840),label='Obsevation')

    selectnum=[]
    for i in range (len(allinall)):
        if allinall[i,9]==14 and allinall[i,8]<=6:
            selectnum.append(i)
    temparray=allinall[selectnum,:]   
    del(selectnum); del(i);
    _,Lx_x_plot,log10_NNsum_plt=plot_Lx_NNsum_gedian(temparray,beaming)
    plt.plot(Lx_x_plot,log10_NNsum_plt,linestyle='--',color=(0.9010,0.550,0.2),label='BH-H')

    selectnum=[]
    for i in range (len(allinall)):
        if allinall[i,9]==14 and allinall[i,8]<=9 and allinall[i,8]>=7:
            selectnum.append(i)
    temparray=allinall[selectnum,:]   
    del(selectnum); del(i);
    _,Lx_x_plot,log10_NNsum_plt=plot_Lx_NNsum_gedian(temparray,beaming)
    plt.plot(Lx_x_plot,log10_NNsum_plt,linestyle=':',color=(0.9010,0.550,0.2),label='BH-He')

    plt.plot(Lx_x, np.log10(N_NSwind + N_NSrb ),linestyle='--',color=(0.4660,0.6740,0.1880),label='NS')

    plt.xlabel('log10(Luminosity) / log10(erg/s)',font1)
    plt.ylabel('log10(Number of ULX)',font1)
    plt.legend(prop=font1)
    plt.yscale('linear')
    plt.xlim((39,40.5))
    plt.ylim((0.5,2))
    #opt.FileName='LuminosityFunctionH-He.jpg';
    #plt.show()


    ##### mass _tb 
    xge = 200
    yge = 100
    temparray = allinall.copy()
    plot_mass_tb_gedian_paper(temparray,xge,yge,' M2-Porb-total');
    #SETplot(' M2-Porb-total')


    temparray=all_rb_1[all_rb_1[i,9]==14,:].copy()
    plot_mass_tb_gedian_paper(temparray,xge,yge,' M2-Porb-BH-RLOF');
    #SETplot(' M2-Porb-BH-RLOF')


    temparray=all_wind_1[all_wind_1[i,9]==14,:].copy()
    plot_mass_tb_gedian_paper(temparray,xge,yge,' M2-Porb-BH-wind')
    #SETplot(' M2-Porb-BH-wind')

    temparray=all_rb_1[all_rb_1[i,9]==13,:].copy() 
    plot_mass_tb_gedian_paper(temparray,xge,yge,' M2-Porb-NS-RLOF')
    #SETplot(' M2-Porb-NS-RLOF')

    temparray=all_wind_1[all_wind_1[i,9]==13,:].copy() 
    plot_mass_tb_gedian_paper(temparray,xge,yge,' M2-Porb-NS-wind')
    #SETplot(' M2-Porb-NS-wind')
    #plt.show()

    ##### plot Lx_porb again
    temparray = allinall.copy() 
    plot_lx_tb_gedian_paper(temparray,xge,yge,' Lx-Porb-total');
    #SETplot(' Lx-Porb-total')


    temparray=all_rb_1[all_rb_1[i,9]==14,:].copy() 
    plot_lx_tb_gedian_paper(temparray,xge,yge,' Lx-Porb-BH-RLOF')
    #SETplot(' Lx-Porb-BH-RLOF')

    
    temparray=all_wind_1[all_wind_1[i,9]==14,:].copy()  

    plot_lx_tb_gedian_paper(temparray,xge,yge,' Lx-Porb-BH-wind')
    #SETplot(' Lx-Porb-BH-wind')


    temparray=all_rb_1[all_rb_1[i,9]==13,:].copy() 

    plot_lx_tb_gedian_paper(temparray,xge,yge,' Lx-Porb-NS-RLOF')
    #SETplot(' Lx-Porb-NS-RLOF')

    selectnum=[]
    for i in range (len(all_wind_1)):
        if all_wind_1[i,9]==13:
            selectnum.append(i)
    temparray=all_wind_1[selectnum,:].copy() 
    del(selectnum); del(i);
    plot_lx_tb_gedian_paper(temparray,xge,100,' Lx-Porb-NS-wind')
    #SETplot(' Lx-Porb-NS-wind')


    plt.show() #记得注释掉！！！！
if __name__ == '__main__':
    PlotAll()