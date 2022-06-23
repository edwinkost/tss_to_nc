#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

#
# PCR-GLOBWB (PCRaster Global Water Balance) Global Hydrological Model
#
# Copyright (C) 2016, Edwin H. Sutanudjaja, Rens van Beek, Niko Wanders, Yoshihide Wada, 
# Joyce H. C. Bosmans, Niels Drost, Ruud J. van der Ent, Inge E. M. de Graaf, Jannis M. Hoch, 
# Kor de Jong, Derek Karssenberg, Patricia López López, Stefanie Peßenteiner, Oliver Schmitz, 
# Menno W. Straatsma, Ekkamol Vannametee, Dominik Wisser, and Marc F. P. Bierkens
# Faculty of Geosciences, Utrecht University, Utrecht, The Netherlands
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import datetime
import time
import re
import glob
import subprocess
import netCDF4 as nc
import numpy as np


def create_nc(input_csv_file = "/scratch/6574882/reanalysis/reanalysis_discharge/pcr_rf_reanalysis_monthly_30arcmin_13.csv", output_nc_file = "test_13.nc" , inp_lon =-31.75, inp_lat = 83.75):
    
    longitude = np.arange(-180.0 + 0.5/2., 180.0,  0.5)
    latitude  = np.arange(  90.0 - 0.5/2., -90.0, -0.5)
    
    # ~ print(long)
    # ~ print(lati)
    
    # ~ long = [inp_lon]
    # ~ lati = [inp_lat]

    nc_format = 'NETCDF4'
    nc_zlib = False
    
    attributeDictionary = {}
    attributeDictionary['institution'] = "Utrecht Univ."
    attributeDictionary['title'      ] = " "
    attributeDictionary['description'] = " "


    # creating nc file

    rootgrp = nc.Dataset(output_nc_file, 'w', format = nc_format)

    #-create dimensions - time is unlimited, others are fixed
    rootgrp.createDimension('lat', len(latitude))
    rootgrp.createDimension('lon', len(longitude))
    rootgrp.createDimension('time', None)

    date_time = rootgrp.createVariable('time', 'f4', ('time',))
    date_time.standard_name = 'time'
    date_time.long_name = 'Days since 1901-01-01'

    date_time.units    = 'days since 1901-01-01'
    date_time.calendar = 'standard'

    lat= rootgrp.createVariable('lat', 'f4', ('lat',))
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'
    lat.standard_name = 'latitude'

    lon= rootgrp.createVariable('lon', 'f4', ('lon',))
    lon.standard_name = 'longitude'
    lon.long_name = 'longitude'
    lon.units = 'degrees_east'

    lat[:] = latitude
    lon[:] = longitude

    # variable name
    varName = "discharge"
    shortVarName = varName
    longVarName  = varName
    standardVarName = varName

    # variable unit
    varUnits = "m3/s"
    
    # missing values
    MV = 1e20
    
    var = rootgrp.createVariable(shortVarName, 'f4', ('time','lat','lon',), fill_value = MV, zlib = nc_zlib)
    var.standard_name = standardVarName
    var.long_name = longVarName
    var.units = varUnits

    for k, v in list(attributeDictionary.items()): setattr(rootgrp,k,v)

    rootgrp.sync()
    rootgrp.close()


    # read the cav file
    with open(input_csv_file) as file_name:
        array = np.loadtxt(file_name, delimiter=",", skiprows = 1, dtype = "str")
    date_string = array[:,0]
    value       = array[:,1].astype(float)     

    
    # write value to a netcdf file
    for i in range(0, len(date_string)):

        timeStamp_string = date_string[i].split('-')
        timeStamp_date = datetime.date(int(timeStamp_string[0]), int(timeStamp_string[1]), int(timeStamp_string[2]))

        # time stamp for reporting
        timeStamp = datetime.datetime(timeStamp_date.year,\
                                      timeStamp_date.month,\
                                      timeStamp_date.day,\
                                      0)
        print(timeStamp)
        
        rootgrp = nc.Dataset(output_nc_file, 'a')

        date_time = rootgrp.variables['time']
        posCnt = len(date_time)
        date_time[posCnt] = nc.date2num(timeStamp, date_time.units,date_time.calendar)

        var_value = value[i]
        
        lat_ind = int(np.where(latitude  == inp_lat)[0])
        lon_ind = int(np.where(longitude == inp_lon)[0])
        
        rootgrp.variables[shortVarName][posCnt,lat_ind,lon_ind] = var_value

        # ~ rootgrp.variables[shortVarName][posCnt, 13, 13] = var_value

        # ~ rootgrp.variables[shortVarName][posCnt, 18, 18] = var_value

        # ~ rootgrp.variables[shortVarName][posCnt,0,0] = var_value

        rootgrp.sync()
        rootgrp.close()

    
def main():

    create_nc(input_csv_file = "/scratch/6574882/reanalysis/reanalysis_discharge/pcr_rf_reanalysis_monthly_30arcmin_13.csv", output_nc_file = "test_13.nc" , inp_lon =-31.75, inp_lat = 83.75)

    # 18,-29.25,83.75
    create_nc(input_csv_file = "/scratch/6574882/reanalysis/reanalysis_discharge/pcr_rf_reanalysis_monthly_30arcmin_18.csv", output_nc_file = "test_18.nc" , inp_lon =-29.25, inp_lat = 83.75)


        
if __name__ == '__main__':
    sys.exit(main())

