#!/share/apps/python/anaconda2.7/bin/python

# Description: This script grabs the basic information about each AR event.
#
# Input: data/ref_data/ARanalysis_ref_data.nc;
#        'data/ARTMIP/AR_identify/%s/%s.%d.nc' % (method, method, year)
#        'data/ARTMIP/MERRA_IVT/IVT.%d.nc' % (year)
#
# Output: 'data/AR_features/part1/%s/ARfeatures_part1.%s.%d.nc' % (method, method, year)
#
# After this script, please aggregate all the files into a single file
#       i.e.,'data/AR_features/part1/AR_feature.part1.%s.1981-2015.nc' % (method)

import numpy as np
import netCDF4 as nc
from skimage.draw import polygon
from skimage import measure
import sys

rootdir = '/raid1/chen423/serdp/archive/GRL2018/'


def get_nc_data(infile, var, proj_option=1):
    #tmpgroup = nc.Dataset(infile, 'r', format='NETCDF4')
    tmpgroup = nc.Dataset(infile, 'r')
    outdata_raw = tmpgroup.variables[var][:]
    tmpgroup.close()
    
    if proj_option==1:
        outdata = np.zeros(outdata_raw.shape)
        outdata[:,:,0:288] = outdata_raw[:,:,288:576]
        outdata[:,:,288:576] = outdata_raw[:,:,0:288]
    if proj_option==2:
        outdata = outdata_raw
    
    return outdata


def int_matrix(inmatrix):
    # since OpenCV only accepts ndarray of uint8 as "figure"
    nx,ny = inmatrix.shape[0:2]
    outdata = np.ones((nx,ny), dtype='uint8')
    for x in np.arange(nx):
        for y in np.arange(ny):
            outdata[x,y] = int(inmatrix[x,y])
            
    return outdata


def crt_AR_event_contours(AR_tag):
    # accepts 2-D arrays, so no time information.
    #AR_tag = int_matrix(AR_sig_data[0,:,:])
    nx, ny = AR_tag.shape[0:2]

    # compute contours
    contours = measure.find_contours(AR_tag, 0.9)

    # number of ARs
    nAR = len(contours)
    
    return nAR, contours


def compute_AR_features(IVT, AR_contour, gridarea, landmask_NA):
    # all the input are 2-D arrays.
    # this mask has 1 for non-AR, 1 for AR, so used for masked_array
    tmp_mask0 = np.ones(IVT.shape)
    #rr, cc = polygon(AR_contour[:,0], AR_contour[:,1], IVT.shape)
    rr, cc = polygon(np.round_(AR_contour[:,0]), np.round(AR_contour[:,1]), IVT.shape)
    tmp_mask0[rr,cc] = 0
    
    # AR full intensity
    tmp_IVT_masked = np.ma.masked_array(IVT, mask=tmp_mask0)
    tmp_ARfull_int = np.ma.average(tmp_IVT_masked, weights=gridarea)
    
    # AR full area
    tmp_area_masked = np.ma.masked_array(gridarea, mask=tmp_mask0)
    tmp_ARfull_area = tmp_area_masked.sum()
    
    # AR land intensity
    mask_AR_NAland = np.maximum(tmp_mask0, landmask_NA)
    tmp_IVT_masked_land = np.ma.masked_array(IVT, mask=mask_AR_NAland)
    tmp_ARland_int = np.ma.average(tmp_IVT_masked_land, weights=gridarea)
    
    # AR land area
    tmp_landarea_masked = np.ma.masked_array(gridarea, mask=mask_AR_NAland)
    tmp_ARland_area = tmp_landarea_masked.sum()
    
    return tmp_ARfull_int, tmp_ARfull_area, tmp_ARland_int, tmp_ARland_area, tmp_mask0


def crt_nc_results(method, year, ARint_full, ARint_coast, ARint_land, ARarea_full, ARarea_land):
    outfile = rootdir+'data/AR_features/part1/%s/ARfeatures_part1.%s.%d.nc' % (method, method, year)
    outgroup = nc.Dataset(outfile, 'w', format='NETCDF4')
    
    # create dimensions
    sdim = outgroup.createDimension('space', 91)
    tdim = outgroup.createDimension('time', None)
    
    # create dimension variables
    svar = outgroup.createVariable('space', 'i4', ('space'))
    svar.long_name = 'Points along the coast'
    svar[:] = np.arange(91)
    tvar = outgroup.createVariable('time', 'i4', ('time'))
    tvar.standard_name = 'time'
    tvar.units = 'hours since %d-01-01 00:00:00' % (year)
    tvar[:] = np.arange(ARint_full.shape[0])*3

    # create variables
    latvar = outgroup.createVariable('lat', 'f4', ('space'))
    latvar.long_name = 'latitude'
    latvar.units = 'degrees_north'
    latvar[:] = west_coast_pt_lats
    
    lonvar = outgroup.createVariable('lon', 'f4', ('space'))
    lonvar.long_name = 'longitude'
    lonvar.units = 'degrees_east'
    lonvar[:] = west_coast_pt_lons
    
    intvar_full = outgroup.createVariable('ARint_full', 'f8', ('time', 'space'))
    intvar_full.long_name = 'mean intensity of AR reaching west coast'
    intvar_full.units = 'kg/m/s'
    intvar_full.note = 'averaged using grid area'
    intvar_full[:] = ARint_full
    
    intvar_coast = outgroup.createVariable('ARint_coast', 'f8', ('time', 'space'))
    intvar_coast.long_name = 'coastal intensity of AR reaching west coast'
    intvar_coast.units = 'kg/m/s'
    intvar_coast[:] = ARint_coast
    
    intvar_land = outgroup.createVariable('ARint_land', 'f8', ('time', 'space'))
    intvar_land.long_name = 'land average intensity of AR reaching west coast'
    intvar_land.units = 'kg/m/s'
    intvar_land.note = 'averaged using grid area'
    intvar_land[:] = ARint_land
    
    areavar_full = outgroup.createVariable('ARarea_full', 'f8', ('time', 'space'))
    areavar_full.long_name = 'area of AR reaching west coast'
    areavar_full.units = 'km^2'
    areavar_full[:] = ARarea_full
    
    areavar_land = outgroup.createVariable('ARarea_land', 'f8', ('time', 'space'))
    areavar_land.long_name = 'land part area of AR reaching west coast'
    areavar_land.units = 'km^2'
    areavar_land[:] = ARarea_land

    outgroup.close()



# get all the reference data
reffile = rootdir+'data/ref_data/ARanalysis_ref_data.nc'
refgroup = nc.Dataset(reffile, 'r', format='NETCDF4')
landmask_NA = refgroup.variables['landmask_NA'][:]
gridarea = refgroup.variables['gridarea'][:]
west_coast_pt_xs = refgroup.variables['pt_x'][:]
west_coast_pt_ys = refgroup.variables['pt_y'][:]
west_coast_pt_lats = refgroup.variables['pt_lat'][:]
west_coast_pt_lons = refgroup.variables['pt_lon'][:]
refgroup.close()



### main function
method  = sys.argv[1]
year = int(sys.argv[2])
proj_option = int(sys.argv[3])
AR_sig_file = rootdir+'data/ARTMIP/AR_identify/%s/%s.%d.nc' % (method, method, year)
AR_sig_data = get_nc_data(AR_sig_file, 'ar_binary_tag', proj_option=proj_option)


IVT_file = rootdir+'data/ARTMIP/MERRA_IVT/IVT.%d.nc' % (year)
IVT_data = get_nc_data(IVT_file, 'IVT')

print(IVT_data.shape)


outdata_meanint = np.zeros((IVT_data.shape[0], 91))
outdata_coastint = np.zeros((IVT_data.shape[0], 91))
outdata_area = np.zeros((IVT_data.shape[0], 91))
outdata_landarea = np.zeros((IVT_data.shape[0], 91))
outdata_landint = np.zeros((IVT_data.shape[0], 91))

for i in np.arange(IVT_data.shape[0]):
    print(i)
    nAR, contours_AR = crt_AR_event_contours(AR_sig_data[i,:,:])
    for n in np.arange(nAR):
        tmp_ARfull_intensity, tmp_ARfull_area, tmp_ARland_intensity, tmp_ARland_area, tmp_mask = compute_AR_features(IVT_data[i,:,:], contours_AR[n], gridarea, landmask_NA)
        for k in np.arange(91):
            if tmp_mask[west_coast_pt_xs[k], west_coast_pt_ys[k]]==0:
                if tmp_ARland_area>0:
                    # AR reaching west coast, now compute the needed quantities
                    # 1. mean intensity over full AR region. Same for each grid
                    outdata_meanint[i,k] = tmp_ARfull_intensity
                    # 2. AR intensity at coast, different for each grid
                    outdata_coastint[i,k] = IVT_data[i,west_coast_pt_xs[k], west_coast_pt_ys[k]]
                    # 3. AR area, full AR region. Same for each grid
                    outdata_area[i,k] = tmp_ARfull_area
                    # 4. AR area, land region. Same for each grid
                    outdata_landarea[i,k] = tmp_ARland_area
                    # 5. AR intensity over land. Same for each grid
                    outdata_landint[i,k] = tmp_ARland_intensity

crt_nc_results(method, year, outdata_meanint, outdata_coastint, outdata_landint, outdata_area, outdata_landarea)
