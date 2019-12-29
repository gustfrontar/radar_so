#!/home/qcradar/.conda/envs/da/bin/python
import datetime as dt                #Datetime operations
import numpy as np                   #Numpy
import os

#=========================================================================================================
# CONFIGURATION SECTION
#=========================================================================================================

#General section
qc_path = "/home/jruiz/share/RADAR/"

datapath = '/home/jruiz/share/DATA/DATOS_RADAR/RMA1/QC_NEW_20181110_RMA1/RMA1/20181110/'  #Main data path.
datapath_out = '/home/jruiz/share/DATA/OBS/LETKF_WRF_SO2KM_RMA1_20181110/'           #Out data path
#deltat = dt.timedelta( seconds=600 )   #Time window (seconds)
#time_offset = 0.0                      #Time offset (from current time)
instrument_list = ['RMA1']  #Instrument list.

file_type_list = ['.H5','.vol','.nc']

remove_local_pkl = False               #Remove intermediate gridded data in pkl format.
remove_local_dat = False               #Remove gridded data in letkf format.
remove_remote_dat = True               #Remove remote letkf files.

#Superobbing section 
output_freq = 60
#        dx    dz   zmax  rmax
grid = [2000, 250,  20e3, 240e3]
opts = {'CZH': [4001, 5, 0], 'CVRAD': [4002, 2]}
#opts = {'CZH': [4001, 5, 0]}
outputpath = datapath_out

c_ini_date='20181110120000'
c_end_date='20181111000000'

#=========================================================================================================
# END OF CONFIGURATION SECTION
#=========================================================================================================

import sys
sys.path.append( qc_path + '/radar_so/')
sys.path.append( qc_path + '/radar_qc/src/python/')

import operational_tools as ot       #Operational tools.
import radar_so as so                #Superobbing module

#Set the dates that will be processed. 
current_date = dt.datetime.utcnow()

#ref_date=dt.datetime(current_date.year, current_date.month, 1, 0, 0, 0)
#freqtimes = deltat.total_seconds()*np.floor((current_date-ref_date).total_seconds()/deltat.total_seconds())
#c_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat   ).strftime('%Y%m%d%H%M%S')
#c_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat*2 ).strftime('%Y%m%d%H%M%S')


print('')
print('=============================================================================')
print('We will process all the files within the following dates:' )
print( c_ini_date )
print( c_end_date )
print('=============================================================================')
print('')

print('')
print('=============================================================================')
print(' GETTING FILE LIST ')
print('=============================================================================')
print('')


#Obtenemos la lista de archivos.
file_list = ot.get_file_list( datapath , c_ini_date , c_end_date , time_search_type='filename' , file_type_list = file_type_list )

print(file_list)

print('')
print('=============================================================================')
print(' READING FILE LIST ')
print('=============================================================================')
print('')

#Obtenemos la lista de objetos radares.
#radar_list = ot.read_multiple_files(  file_list , instrument_list )

for my_file_radar in file_list :

      print('')
      print('=============================================================================')
      print(' SUPEROBBING')
      print('=============================================================================')
      print('')

      #Call SO routine 
      radar = ot.read_multiple_files(  [my_file_radar] , instrument_list )[0]

      letkf_filelist = so.main_radar_so(radar, output_freq, grid, opts, datapath_out  )

      for my_file in letkf_filelist :

          my_time_datetime = ot.get_time_from_filename( my_file )

          my_time = dt.datetime.strftime( my_time_datetime  + dt.timedelta( seconds=3600) , '%Y%m%d_%H' )

          my_minute = dt.datetime.strftime( my_time_datetime , '%M' )
 
          complete_path = datapath_out + '/radar/' + my_time

          if not os.path.isdir( complete_path )  :

             os.makedirs( complete_path )

          os.system('ln -sf ' + my_file + ' ' + complete_path + '/' + os.path.basename(my_file) )

          if my_minute == '00'   :  #Copy the file in the next hour folder as well.

             my_time = dt.datetime.strftime( my_time_datetime , '%Y%m%d_%H' )

             complete_path = datapath_out + '/radar/' + my_time

             if not os.path.isdir( complete_path )  :

                os.makedirs( complete_path )

             os.system('ln -sf ' + my_file + ' ' + complete_path + '/' + os.path.basename(my_file) )
             




