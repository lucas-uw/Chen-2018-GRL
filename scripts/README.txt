Note on the script:

1. Within each script, there is a variable "rootdir". This is where you download this whole repo and put on your local machine.

2. Under the rootdir, put the yearly split AR tag files at 'data/ARTMIP/AR_identify/%s/%s.%d.nc' % (method, method, year)

3. Under the rootdir, put the MERRA2 IVT files at 'data/ARTMIP/MERRA_IVT/IVT.%d.nc' % (year)

4. after step01, please aggregate all the files into a single file
       i.e.,'data/AR_features/part1/AR_feature.part1.%s.1981-2015.nc' % (method)
       



About the figures each script creates:
step03.  Figure S7, S8
step04.
step05.  Figure 4, 5
step06.  Figure 2, 3
step07.  Figure S5, S6
step08.  Figure S9, S10, S11
step09.  Figure S2
step10.  Figure S3, S4
step11.  Figure 1
