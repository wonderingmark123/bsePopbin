# bsePopbin

There are 4 file folders in this project.

## Fortran0415

modified fortran code BSE
- popbin_mine.f: main function for POPBIN program, get the evolution details of every stage.
- popbin_hanchen.f: POPBIN program. This program is for calculating the birth rate of each binary systems when the Neutron star born.
More details are available in README_BSE, README_SSE and README_NEW.
The history of modified could be found in [note.txt](./Fortran0415/note.txt) here.

## FortranOri

This is the original BSE code for backup.

## Python

This file folder contains all the analysis python program. More details could be found in its [README.md](./python/README.md)

## Data

This is the file folder for saving data and will not be uploaded to github.

### PaperPlotData

This is the file folder for saving mesa data and popbin data with label being ``1117``.

- 1117_paper.mat: popbin data with label being ``1117``
- hanchen_first_neutron1124.out: popbin result for mesa plot.
- Lx_N1207.npy: selected data from mesa output data.
  
### Contributor

    Haotian Song modified in 2021/5/7
    Email: 490785554@qq.com 



