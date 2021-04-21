# Binary Population Synthesis Code

## Usage

All the functions are available in [utils.py](utils.py) which are not necessary to be modified for general plots.

[main.py](main.py) is the main function and some examples for plots. And make sure that use .copy() to copy the numpy array during the selection. Otherwise, the original numpy array will be modified in some functions.

### [main.py](main.py)
- main0415: Loading data from Fortran Program Output with version later than 0415. More details are available in [Fortran0415/popbin_mine.f](../Fortran0415/popbin_mine.f)
- main1117: Loading data from previous versions, also an example for loading data from MATLAB data files(.mat)
- [PlotAll](paper_plot.py): Plot all the figures in this function. You can also plot individual figure in main function. More details are available in example. 
- ImportOriData: 
### Contributor

> Haotian Song modified in 2021/4/21
> Email: 490785554@qq.com