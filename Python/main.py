from utils import *
import numpy as np
import os
from paper_plot import massFunction,PlotAll

def main0415():
    print('Loading data now')
    all_rb_1,all_wind_1 = ImportOriData(
    NumParallel     = 8,
    DataFolder      = '..\Data',
    OutName         = '0415',
    SaveNPY         = False,
    LoadNPY         = True)
    PlotAll(all_rb_1,all_wind_1,
    FileFolder = 'D:/study/bsegit/bsePopbin/Data/PaperPlotData')
    
    return

def main1117():
    PlotAll([],[],
    FileFolder = 'D:/study/bsegit/bsePopbin/Data/PaperPlotData')
    return
if __name__ == '__main__':
    main0415()