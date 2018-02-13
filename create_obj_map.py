#
import sys, os, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as P
from optparse import OptionParser
import datetime as dtime
import xarray as xr
import matplotlib.cm as cm
from datetime import datetime
import datetime as DT
import netCDF4 as ncdf

import pyresample
from scipy.ndimage import gaussian_filter
from skimage.morphology import watershed
from skimage.feature import peak_local_max
from skimage.measure import regionprops
from scipy import ndimage as ndi

# Parameters

_innov_min    = 35.
_level        = 10
_max_obj_size = (10,10)
_obj_inflate  = 3

time_format = "%Y%m%d%H%M"
dart_calendar = datetime.strptime("160101010000", "%Y%m%d%H%M")

#-------------------------------------------------------------------------------
#
def obs_seq_read_netcdf(filename, retFileAttr = False):
    if retFileAttr == False:
        return xr.open_dataset(filename).to_dataframe()
    else:
        xa = xr.open_dataset(filename)
        return xa.to_dataframe(), xa.attrs

#-------------------------------------------------------------------------------
#
def obs_seq_get_obtype(df, kind=None, name=None):

    if kind:
        idx = df['kind'] == kind
        return df[idx]

    if name:
        idx = df['name'] == name
        return df[idx]

    print("\n OBS_SEQ_GET_OBTYPE:  no kind or name specified, exiting \n")
    sys.exit(-1)

#=========================================================================================
# Using innovations from posterior obs_seq files, create a list of objects where
#     the innov > _innov_min
#-------------------------------------------------------------------------------
# Main function defined to return correct sys.exit() calls

def main(argv=None):

    if argv is None:
           argv = sys.argv

# Command line interface for DART_cc

    parser = OptionParser()

    parser.add_option("-f", "--file",  dest="file",  default=None, type="string",
                  help = "Obs_seq_final file to convert")

    parser.add_option("-w", "--wrf",  dest="wrf_file",  default=None, type="string",
                  help = "WRF file to specify the domain, any one of them will do")

    parser.add_option("-i", "--interactive", dest="interactive", default=False,  action="store_true",     \
                     help = "Boolean flag to dump infomation and plot image to screen")

    (options, args) = parser.parse_args()

    if options.wrf_file == None or options.file == None:

        print("\n CREATE_OBJ_MAP:  You did not specify either the WRF domain file or the netCDF obs_seq file!!" )
        print("  .............EXITING.... \n ")
        parser.print_help()
        print
        sys.exit(1)

    else:
        filenames = [options.file, options.wrf_file]
        date = datetime.strptime(options.file[-15:-3], "%Y%m%d%H%M")
        day  = datetime.strptime(options.file[-15:-7], "%Y%m%d")
        days = (date - dart_calendar).days
        secs = (date - day).total_seconds()

        output_file = ("refl_obs_%d_%d.txt" % (np.floor(secs), days))
        print("\n CREATE_OBJ_MAP: writing to txt file:  %s" % output_file)

    dataset, fileAttrs = obs_seq_read_netcdf(filenames[0], retFileAttr = True)
    kind = fileAttrs['RADAR_REFLECTIVITY']
    field = obs_seq_get_obtype(dataset, kind=kind)

    index = field.innov > _innov_min
    
    if len(index) > 0:
        innov    = field[index].innov.values
        lats    = field[index].lat.values
        lons    = field[index].lon.values

# Exit if there are no objects

    if len(innov) < 1:
        print("\n\n  CREATE_OBJ_MAP:  No innovations found > %d, exiting\n\n" % _innov_min)
        sys.exit(-1)

    wrf_grid_file = ncdf.Dataset(filenames[1])

    glats = wrf_grid_file['XLAT'][0]
    glons = wrf_grid_file['XLONG'][0]

    nx = glats.shape[1]
    ny = glats.shape[0]

    innov_2d = np.zeros(glats.shape)

    grid  = pyresample.geometry.GridDefinition(lats=glats, lons=glons)
    swath = pyresample.geometry.SwathDefinition(lons=lons, lats=lats)

# Determine nearest (w.r.t. great circle distance) neighbour in the grid.
    _, _, index_array, distance_array = pyresample.kd_tree.get_neighbour_info(
        source_geo_def=grid, target_geo_def=swath, radius_of_influence=50000,
        neighbours=1)

#get_neighbour_info() returns indices in the flattened lat/lon grid. Compute the 2D grid indices:

    index_array_2d = np.unravel_index(index_array, grid.shape)

    for n in np.arange(len(index_array)):
        innov_2d[index_array_2d[0][n], index_array_2d[1][n]] = innov[n]

    innov_2d = gaussian_filter(innov_2d, sigma=1.0)

    innov_binary = np.where(innov_2d > 1., 1, 0)

    distance = ndi.distance_transform_edt(innov_binary)
    local_max = peak_local_max(distance, indices=False,  footprint=np.ones(_max_obj_size), labels=innov_binary)
    markers = ndi.label(local_max)[0]
    var_labels = watershed(-distance, markers, mask=innov_binary)
    var_props = regionprops(var_labels, innov_2d)

    outfile = open(output_file, "w")

# To make things more "fuzzy", we are going to then re-inflate the physical size of the
# objects to be "obj_size" by creating a 3D cube around each object where noise will be
# added.  This still will reduce the number of noise points by an ~ order of magnitude
# Hardcoding this now for simplicity to a 3x3x3 cube

    for k in range(1+len(var_props)-1):
        
        ii = int(np.round(var_props[k].centroid[1]))
        jj = int(np.round(var_props[k].centroid[0]))
        kk = _level

        iip1 = min(ii+1, nx)
        iim1 = max(ii-1, 1)
        jjp1 = min(jj+1, ny)
        jjm1 = max(jj-1, 1)
        kkm1 = max(kk-1, 1)
        kkp1 = min(kk+1, 50)
        
        if options.interactive:
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (ii,   jj,   kkm1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iip1, jj,   kkm1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iim1, jj,   kkm1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (ii,   jjp1, kkm1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (ii,   jjm1, kkm1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iip1, jjp1, kkm1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iim1, jjp1, kkm1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iip1, jjm1, kkm1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iim1, jjm1, kkm1, innov_2d[jj,ii]))

            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (ii,   jj,   kk,   innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iip1, jj,   kk,   innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iim1, jj,   kk,   innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (ii,   jjp1, kk,   innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (ii,   jjm1, kk,   innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iip1, jjp1, kk,   innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iim1, jjp1, kk,   innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iip1, jjm1, kk,   innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iim1, jjm1, kk,   innov_2d[jj,ii]))
            
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (ii,   jj,   kkp1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iip1, jj,   kkp1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iim1, jj,   kkp1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (ii,   jjp1, kkp1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (ii,   jjm1, kkp1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iip1, jjp1, kkp1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iim1, jjp1, kkp1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iip1, jjm1, kkp1, innov_2d[jj,ii]))
            print("   %3.3i    %3.3i    %3.3i    %4.2f\n" % (iim1, jjm1, kkp1, innov_2d[jj,ii]))

        outfile.write("   %i    %i    %i    %4.2f\n" % (ii,   jj,   kkm1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iip1, jj,   kkm1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iim1, jj,   kkm1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (ii,   jjp1, kkm1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (ii,   jjm1, kkm1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iip1, jjp1, kkm1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iim1, jjp1, kkm1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iip1, jjm1, kkm1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iim1, jjm1, kkm1, innov_2d[jj,ii]))

        outfile.write("   %i    %i    %i    %4.2f\n" % (ii,   jj,   kk,   innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iip1, jj,   kk,   innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iim1, jj,   kk,   innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (ii,   jjp1, kk,   innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (ii,   jjm1, kk,   innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iip1, jjp1, kk,   innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iim1, jjp1, kk,   innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iip1, jjm1, kk,   innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iim1, jjm1, kk,   innov_2d[jj,ii]))
            
        outfile.write("   %i    %i    %i    %4.2f\n" % (ii,   jj,   kkp1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iip1, jj,   kkp1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iim1, jj,   kkp1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (ii,   jjp1, kkp1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (ii,   jjm1, kkp1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iip1, jjp1, kkp1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iim1, jjp1, kkp1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iip1, jjm1, kkp1, innov_2d[jj,ii]))
        outfile.write("   %i    %i    %i    %4.2f\n" % (iim1, jjm1, kkp1, innov_2d[jj,ii]))

    outfile.close()

# Plot out the innovations mapped to grid, and objects

    fig, axes = P.subplots(1, 2, sharey=True, figsize=(18,8))
    axes[0].contourf(innov_2d,cmap=P.cm.spectral, interpolation='nearest')
    axes[0].set_title("Innovation Map: %s" % date.strftime("%Y-%m-%d %H:%M"))
    #axes[0].set_aspect("equal")
    
    axes[1].contourf(var_labels, cmap=P.cm.spectral, interpolation='nearest')
    axes[1].set_title("Object Map: %s" % date.strftime("%Y-%m-%d %H:%M"))
    #axes[1].set_aspect("equal")
    for k in range(1+len(var_props)-1):
        axes[1].plot(int(np.round(var_props[k].centroid[1])), int(np.round(var_props[k].centroid[0])), "o", markersize=2.0, color='w')

    P.savefig(output_file[:-4]+".png")
    if options.interactive:
        P.show()


#-------------------------------------------------------------------------------
# Main program for testing...
#
if __name__ == "__main__":
    sys.exit(main())

# End of file
