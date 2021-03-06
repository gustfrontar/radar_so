"""
Functions to perform radar superobbing for LETKF-WRF data assimilation system
"""
from collections import defaultdict
from datetime import datetime, timedelta
import pickle
import struct
import numpy as np
import pyart
import matplotlib.pyplot as plt
import sys 
import os
sys.path.append( '/home/jruiz/share/RADAR/radar_so/fortran/' )
import ntpath
from cs  import cs   #Fortran code routines.

class SoFields(object):
    """
    Superobbed radar fields.

    Parameters
    ----------
    radar : Radar
        PyArt-Like radar object
    dx, dz : int
        Cartesian grid horizontal and vertical resolution, respectively, in
        meters.
    max_height: int
        Cartesian grid maximum height in meters.
    mindbz : int
        Reflectivity minimum value.

    Attributes
    ----------
    radar : Radar
        PyArt-Like radar object
    grid : SoGrid
        Cartesian grid object
    fields : dict of dicts
        Superobbed fields (i.e. reflectivity and Doppler velocity), containing
        the field's data, range, azimuth and elevation box-averaged values and
        the number of data points used to perform the average.

    Methods
    -------
    """
    def __init__(self, input, variables, rays, grid, opts):
        print('SO Fields object')
        self._options = opts
        self.radar = SoRadar(input, variables, rays)
        self.grid = SoGrid(self.radar, grid[0], grid[1], grid[2],grid[3])
        self.fields = defaultdict(dict)

        self.compute_so(variables)

    def compute_so(self, variables):
        """
        Compute superobbing for certain variables

        Parameters
        ----------
        variables : list of str
            Radar variable names as in radar file.
        """
        vars_data = {}
        for var in variables:
            if self._check_field_exists(var):
                # Allocate vars in fields attribute
                self.allocate_var(var, self.grid.nlon, self.grid.nlat, self.grid.nlev)
                # Get data from radar object
                vars_data[var] = self.radar.fields[var]

                # Convert dBZ to power
                #if 'reflectivity' in vars_data[var]['standard_name']:
                if var == 'CZH':
                    #tmp_mask = np.ones( np.shape( vars_data[var]['data'] ) ).astype(bool)
                    #tmp_mask[ vars_data[var]['data'] == vars_data[var]['_FillValue']  ] = False
                    vars_data[var]['data'] = np.power(10.0, vars_data[var]['data']/10.0)
                    vars_data[var]['data'].data[ vars_data[var]['data'].mask  ] = vars_data[var]['_FillValue']


            #print(vars_data[var]['data'].shape)

        # Average data
        self.compute_grid_boxmean(vars_data)

    def allocate_var(self, var, nx, ny, nz):
        """ Allocate variables to a key in self.fields """
        self.fields['grid_' + var]['nobs'] = np.zeros((nz, ny, nx))
        self.fields['grid_' + var]['data'] = np.zeros((nz, ny, nx))
        self.fields['grid_' + var]['az'] = np.zeros((nz, ny, nx))
        self.fields['grid_' + var]['el'] = np.zeros((nz, ny, nx))
        self.fields['grid_' + var]['ra'] = np.zeros((nz, ny, nx))
        self.fields['grid_' + var]['id'] =  self._options[var][0]
        self.fields['grid_' + var]['error'] =  self._options[var][1]
        if len(self._options[var]) == 3:
            self.fields['grid_' + var]['min'] = self._options[var][2]

    def compute_grid_boxmean(self, variables):
        import numpy.matlib 
        """
        Compute the i, j, k for each radar grid point and box-average the data.

        Parameters
        ----------
        variables : list of numpy array
        """

        nr=self.radar.ngates
        na=self.radar.nrays
        nvar=len( self.fields.keys() )

        lon_ini=self.grid.lon[0, 0]
        lat_ini=self.grid.lat[0, 0]
        z_ini=0.0

        dlon=self.grid.dlon
        dlat=self.grid.dlat
        dz  =self.grid.dz

        nlon=self.grid.nlon
        nlat=self.grid.nlat
        nz  =self.grid.nlev


        local_undef = -8888.0

        #Group all variables that will be "superobbed" into one single array.
        #This is to take advantage of paralellization. Different columns will be processed in parallel.
        datain  =np.zeros( ( na*nr , 4*nvar ) )
        latin   =np.zeros( ( na*nr  ) )
        lonin   =np.zeros( ( na*nr  ) )
        zin     =np.zeros( ( na*nr  ) )

        var_names = []

        weigth=np.zeros( 4*nvar ).astype(bool)  #This flag decides wether the superobbing will be weigthed by the reflectivity or not (for each variable)
        weigth_ind=np.zeros( 4*nvar )
        weigth_vars=False  #TODO this should go to the configuration. (Enable or disable reflectivity weithning)
        w_i=0


        for iv , var in enumerate( variables )  :

            for ir in range( nr ) :
               #for ir in range( nr ) :
               datain[ir*na:(ir+1)*na,4*iv+0]=self.radar.azimuth['data'][:]
               datain[ir*na:(ir+1)*na,4*iv+1]=self.radar.elevation['data'][:]
               datain[ir*na:(ir+1)*na,4*iv+2]=self.radar.range['data'][ir]
               datain[ir*na:(ir+1)*na,4*iv+3]=variables[var]['data'].data[:,ir]

               latin[ir*na:(ir+1)*na] = self.radar.gate_latitude['data'][:,ir]
               lonin[ir*na:(ir+1)*na] = self.radar.gate_longitude['data'][:,ir]
               zin  [ir*na:(ir+1)*na] = self.radar.gate_altitude['data'][:,ir]

            #Change for variable undef into local_undef.
            tmp_mask = ( datain[:,4*iv+3] == variables[var]['_FillValue'] )
            for ii in range( 0 , 4 ) :
               datain[ tmp_mask , 4*iv+ii ] = local_undef

            if var == 'CHZ' and weigth_vars :   #TODO Replace CHZ by the "reflectivity name"
               w_i = 4*iv+3 #Reflectivity is the variable that can be used as a weight.
               weigth=np.ones( 4*nvar ).astype(bool) 
               weigth[w_i-4:w_i]=False                 #Do not weigth the reflectivity.
               weigth_ind=np.ones( 4*nvar )*w_i

        is_angle=np.zeros( 4*nvar )             #This flag decides if variable will be treated as an angle in the superobbing (usually true for the azimuth)
        is_angle[range(0,4*nvar,4)]=True        #Azimuths will be treated as angles.

        #The function outputs 
        #data_ave - weigthed average
        #data_max - maximum value
        #data_min - minimum value
        #data_std - standard deviation (can be used to filter some super obbs)
        #data_n   - number of samples  (can be used to filter some super obbs)
        #data_w   - sum of weigths

        print('Doing average')
        [data_ave , data_max , data_min , data_std , data_n , data_w ]=cs.com_interp_boxavereg(
                             xini=lon_ini,dx=dlon,nx=nlon,yini=lat_ini,dy=dlat,ny=nlat,
                             zini=z_ini  ,dz=dz  ,nz=nz  ,nvar=4*nvar,
                             xin=lonin,yin=latin,zin=zin,datain=datain,undef=local_undef,nin=na*nr,   
                             weigth=weigth,weigth_ind=weigth_ind,is_angle=is_angle)

        #"Unpack" the superobbed data.
        for iv , var in enumerate( variables )  :
            self.fields['grid_' + var]['az']=data_ave[:,:,:,0+iv*4]
            self.fields['grid_' + var]['el']=data_ave[:,:,:,1+iv*4]
            self.fields['grid_' + var]['ra']=data_ave[:,:,:,2+iv*4]
            self.fields['grid_' + var]['data']=data_ave[:,:,:,3+iv*4]
            self.fields['grid_' + var]['nobs']=data_n[:,:,:,3+iv*4]

        #print('MAX AZ', np.max( self.fields['grid_' + var]['az']  ) )

    #*********************
    # Auxiliary functions
    #*********************
    def _radar2grid(self, ray, range):
        reali = (self.radar.gate_longitude['data'][ray, range]-self.grid.lon[0, 0])/self.grid.dlon
        realj = (self.radar.gate_latitude['data'][ray, range]-self.grid.lat[0, 0])/self.grid.dlat
        realk = (self.radar.gate_altitude['data'][ray, range]-self.grid.z[0, 0, 0])/self.grid.dz

        i = np.round(reali).astype('int')
        j = np.round(realj).astype('int')
        k = np.round(realk).astype('int')

        return [k, j, i]

    def _check_point_in_grid(self, i, j, k):
        """ Check if a gridpoint is inside grid domain """
        return True if 0 <= k < self.grid.nlev and \
            0 <= i < self.grid.nlon and \
            0 <= j < self.grid.nlat else False

    def _check_field_exists(self, field_name):
        if field_name not in self.radar.fields.keys():
            raise KeyError('Field not available: ' + field_name)
        return True

class SoRadar(object):
    """
    A PyArt-Like radar object with necessary attributes and variables.

    Parameters
    ----------
    filename : str
        Radar file.
    variables: list
        Variables to extract from radar object.
    ray_interval : list
        Initial and final ray number to reduce radar fields data.

    Attributes
    ----------
    gate_longitude, gate_latitude : LazyLoadDict
        Geographic location of each gate.  The projection parameter(s) defined
        in the `projection` attribute are used to perform an inverse map
        projection from the Cartesian gate locations relative to the radar
        location to longitudes and latitudes.

    """
    def __init__(self, input, variables, ray_interval):
        self.__radar = input
        self.__variables = variables
        self.__ray_interval = ray_interval

        self.gate_longitude = {}
        self.gate_latitude = {}
        self.gate_altitude = {}
        self.azimuth = {}
        self.elevation = {}
        self.fields = defaultdict(dict)

        self.ngates = self.__radar.ngates
        self.nrays = ray_interval[-1] - ray_interval[0]
        self.longitude = self.__radar.longitude
        self.latitude = self.__radar.latitude
        self.altitude = self.__radar.altitude
        self.range = self.__radar.range

        self.get_data()

    def get_data(self):
        self.gate_longitude['data'] = self._extract_rays(self.__radar.gate_longitude['data']).copy()
        self.gate_latitude['data'] = self._extract_rays(self.__radar.gate_latitude['data']).copy()
        self.gate_altitude['data'] = self._extract_rays(self.__radar.gate_altitude['data']).copy()
        self.azimuth['data'] = self._extract_rays(self.__radar.azimuth['data']).copy()
        self.elevation['data'] = self._extract_rays(self.__radar.elevation['data']).copy()

        for var in self.__variables:
            self.fields[var] =  self.__radar.fields[var].copy()
            self.fields[var]['data'] = self._extract_rays(self.__radar.fields[var]['data']).copy()

    def _extract_rays(self, a):
        if a.ndim == 1:
            return a[self.__ray_interval[0]:self.__ray_interval[-1]]
        return a[self.__ray_interval[0]:self.__ray_interval[-1], :]

class SoGrid(object):
    """
    A cartesian grid.

    Parameters
    ----------
    radar : Radar
        Pyart radar object to use for georeference.
    dx : int
        Horizontal grid dimension in meters.
    dz : int
        Vertical grid dimension in meters.
    maxz : int
        Maximum height in meters.

    Attributes
    ----------
    dx : int
        Horizontal grid dimension in meters.
    dz : int
        Vertical grid dimension in meters.
    nlon, nlat, nlev : int
        Dimensions of the longitude, latitude and height axis, respectively.
    dlon, dlat : float
        Distance in degrees between to grid points in the longitude and latitude
        axis, respectively.
    lon : numpy array
        Longitude for each grid point.
    lat : numpy array
        Latitude for each grid point.
    z : numpy array
        Height for each grid point.
    """
    def __init__(self, radar, dx, dz, maxz , maxrange = None ):

        self.dx = dx
        self.dz = dz

        # Compute possible value for `nlon` in order to cover the maximum radar range
        if maxrange == None  :
            maxrange = np.max(radar.range['data'])
        self.nlon = np.ceil(2.*maxrange/self.dx).astype('int')
        self.nlat = self.nlon
        self.nlev = np.ceil(maxz/self.dz).astype('int')

        # Force grid dimension to be odd
        if np.mod(self.nlon, 2) == 0:
            self.nlon += 1
        if np.mod(self.nlat, 2) == 0:
            self.nlat += 1

        self.dlon = None
        self.dlat = None
        self.lon = np.zeros((self.nlat, self.nlat))*np.nan
        self.lat = np.zeros((self.nlat, self.nlon))*np.nan
        self.z = np.zeros((self.nlev, self.nlat, self.nlon))*np.nan

        radar_lon = radar.longitude['data']
        radar_lat = radar.latitude['data']

        self.dx2ddeg(radar_lat)
        self.lat_lon(radar_lon, radar_lat)
        self.zlevs()

    def dx2ddeg(self, radar_lat):
        """
        Translate `dx` into an appropiate delta in longitude and latitude
        Hopefully, nobody will put a radar at the pole.
        """
        re = 6371e3
        self.dlon = float(np.rad2deg(self.dx/(re*np.cos(np.deg2rad(radar_lat)))))
        self.dlat = float(np.rad2deg(self.dx/re))

    def lat_lon(self, radar_lon, radar_lat):
        """ Define latitudes and longitudes """
        for i in range(self.nlat):
            for j in range(self.nlon):
                self.lon[i, j] = radar_lon + self.dlon*(-1.0-(self.nlon-1.0)/2.0 + (j+1))
                self.lat[i, j] = radar_lat + self.dlat*(-1.0-(self.nlat-1.0)/2.0 + (i+1))

    def zlevs(self):
        """ Define height levels """
        for k in range(self.nlev):
            self.z[k, :, :] = k*self.dz

def main_radar_so(input, output_freq, grid_dims, options ,outputpath=None):
    """
    Perform superobbing to radar volume

    Parameters
    ----------
    filename : str
        Radar file.
    output_freq : int, in seconds
        LETKF files output frequency in seconds.
    grid_dims : list [dx, dz, max_height]
        Cartesian grid horizontal and vertical dimensions and maximum height, in meters.
    radar : objecto opcional.
    options : dict
        Information needed to write the binary file in LETKF format.
        {'ref': [id, error, minref], 'dv': [id, error]}
    """
    # Read radar volume using pyart
    if isinstance(input, str):
        radar = pyart.io.read(filename)
    else:
        radar = input
    #print(radar)

    if outputpath == None  :
        outputpath='./'
    os.makedirs( outputpath + '/grid/' ,exist_ok=True)
    os.makedirs( outputpath + '/letkf/',exist_ok=True)

    # Get reflectivity and Doppler velocity variable name
    vars_name = get_vars_name(radar,options)
    print(vars_name)
    #vars_name = radar.fields.keys()

    # Get radar start and end times
    inidate, enddate = get_dates(radar)
    # Create output files list according to out_freq
    output_dates = get_letkf_outputs(inidate, enddate, output_freq)
    #Loop over output files
    iray = 0
    inirayidx = 0
    outfile_list = []

    if not vars_name  :

       return outfile_list

    for mydate in output_dates :

        # Get radar rays that contribute to current date
        top_second = ( ( mydate + timedelta(seconds=output_freq/2.0) ) - inidate ).total_seconds()
        bot_second = ( ( mydate - timedelta(seconds=output_freq/2.0) ) - inidate ).total_seconds()

        my_rays= np.squeeze( np.where( np.logical_and( radar.time['data'] >= bot_second , radar.time['data'] < top_second ) ) ) 

        #diff = (output_dates[date][1]-inidate).total_seconds()
        #while iray < radar.nrays and np.abs(radar.time['data'][iray] - diff) > 1e-3  :
        #    iray += 1
        #endrayidx = iray
        #ray_limits = [inirayidx, endrayidx]
        if np.size( my_rays ) > 1 :
           ray_limits = [ my_rays[0] , my_rays[-1] ]
           #print( ray_limits )

           # Compute superobbing
           so = SoFields(radar, vars_name, ray_limits, grid_dims, options)

           #print('PLOTTING ORIGINAL')
           #display = pyart.graph.RadarDisplay(radar)
           #for sweep in range(radar.nsweeps):
           #   fig = plt.figure()
           #   display.plot_ppi('CZH', sweep=sweep)
           #   plt.show()

           '''
           print('PLOTTING SO')
           for lev in so.fields['grid_CZH']: #sso.grid.nlev):
              if lev != 'id' and lev != 'error' and lev != 'min':
                  print(lev)
                  fig = plt.figure()
                  plt.pcolormesh(so.fields['grid_CZH'][lev][3,:,:])
                  plt.colorbar()
                  plt.show() 
           '''

           # Check if there is an exisiting file to update the box average
           tmpfile = outputpath + '/grid/' + radar.metadata['instrument_name'] + '_' + date2str(mydate) + '.pkl'
        
           if check_file_exists(tmpfile):
               print('Updating boxmean from previous file ' + tmpfile)
               tmp_so = load_object(tmpfile)
               update_boxaverage(tmp_so, so)

           '''
           print('PLOTTING SO')
           for lev in so.fields['grid_CZH']: #sso.grid.nlev):
              if lev != 'id' and lev != 'error' and lev != 'min':
                  print(lev)
                  fig = plt.figure()
                  plt.pcolormesh(so.fields['grid_CZH'][lev][3,:,:])
                  plt.colorbar()
                  plt.show()
           '''

           # Write intermediate file
           write_object(tmpfile, so.fields)

           # Write LETKF file
           outfile = outputpath + '/letkf/' + radar.metadata['instrument_name'] + '_' + date2str(mydate) + '.dat'
           outfile_list.append( outfile )
           write_letkf(outfile, so)

           inirayidx = iray
     
    return outfile_list 

def get_output_dates_list(dates):
    '''
    dates : dict of datetime objects
    '''
    datelist = []
    for key in dates:
       datelist.append(date2str(key))
    return datelist


def get_vars_name(obj,options):
    """ Get reflectivity and Doppler veolocity names from radar object """
    variable_name = []
    for key in obj.fields.keys()  :
        if key in options.keys()  :
              #if key != 'CM' and key != 'TV' and ('reflectivity' in obj.fields[key]['standard_name'] or \
              #    'radial_velocity' in obj.fields[key]['standard_name']):
           variable_name.append(key)
           #print( variable_name )
    return variable_name

def get_dates(obj):
    """ Get radar initial and final dates from radar object """
    initime = obj.time['units'].split(' ')[2]
    initial_date = str2date(initime)
    delta = timedelta(seconds=obj.time['data'][-1].tolist())
    end_date = initial_date + delta
    return initial_date, end_date

def get_letkf_outputs(inidate, enddate, freq):
    """
    Create a list of letkf output filenames
    initime, endtime : datetime objects
    out_freq : int in seconds
    """
    ref_date =  datetime(inidate.year, inidate.month, 1, 0, 0, 0)

    #if freq < 3600:
    # Get first possible output date
    #inirounddate = datetime(inidate.year, inidate.month, inidate.day, inidate.hour, 0, 0)
    freqtimes = freq*round((inidate-ref_date).total_seconds()/freq)
    outputini = ref_date + timedelta(seconds=freqtimes)

    # Get last possible output date
    
    #endrounddate = datetime(inidate.year, inidate.month, inidate.day, inidate.hour+1, 0, 0)
    freqtimes = freq*round((enddate-ref_date).total_seconds()/freq)
    outputend = ref_date + timedelta(seconds=freqtimes)

    # Create output list
    delta = timedelta(seconds=freq)
    output_dates = defaultdict(list)
    for date in datespan(outputini, outputend+delta, delta):
        iniinterval = date - timedelta(seconds=freq/2)
        endinterval = date + timedelta(seconds=freq/2)
        output_dates[date] = [iniinterval, endinterval]

    # Check radar intial and final dates are in output list
    check_date_in_interval(inidate, output_dates[outputini][0], output_dates[outputini][1])
    check_date_in_interval(enddate, output_dates[outputend][0], output_dates[outputend][1])

    return output_dates
    #else:
    #    raise ValueError('Not a valid output frequency:' + freq)

def check_date_in_interval(date, lower, upper):
    if not lower <= date <= upper:
        raise ValueError('Date not in interval: ' + date)
    return True

def check_file_exists(filename):
    """ file : pickle file """
    import os
    return True if os.path.exists(filename) else False

def update_boxaverage(old_obj, new_obj):
    """ Update superobbing object with previous data """
    for key in new_obj.fields.keys():
        #Check if the current key is present in the old_object.
        #If so merge the two objects.
        if key in old_obj.keys()  :
           nobs_old = old_obj[key]['nobs']
           nobs_new = new_obj.fields[key]['nobs']
           nobs_tot = nobs_old + nobs_new

           for subkey in new_obj.fields[key].keys():
               if subkey != 'nobs' and subkey != 'id' and subkey != 'error' and subkey != 'min':
                   #print( nobs_old.shape , old_obj[key][subkey].shape )
                   new_obj.fields[key][subkey] = \
                       nobs_new*new_obj.fields[key][subkey] + nobs_old*old_obj[key][subkey]
                   new_obj.fields[key][subkey][ nobs_tot > 0 ] = new_obj.fields[key][subkey][ nobs_tot > 0 ] / nobs_tot[ nobs_tot > 0 ]
                   #   np.ma.masked_invalid((nobs_new*new_obj.fields[key][subkey] +\
                   #   nobs_old*old_obj[key][subkey])/np.ma.masked_invalid(nobs_tot))
           new_obj.fields[key]['nobs'] = nobs_tot

def write_object(filename, obj):
    print('WRITING PICKLE FILE ' + filename)
    with open(filename, 'wb') as fileout:  # Overwrites any existing file.
        pickle.dump(obj, fileout, pickle.HIGHEST_PROTOCOL)

def load_object(filename):
    with open(filename, 'rb') as filein:
        return pickle.load(filein)

def write_letkf(filename, obj):
    print('WRITING BINARY LETKF FILE ' + filename)
    tmp1 = np.array([4])
    tmp2 = tmp1*7
    nobs = 0
    wk = np.empty(7)

    nvar = len( obj.fields.keys() )

    tmp_error=np.zeros( nvar )
    tmp_id   =np.zeros( nvar )
    tmp_lambda =  12.0 #TODO check this value and get it from the radar object.
                       #Temporarily this value is assumed to be equal to S-BAND 
                       #an specific C-BAND operator has not been coded yet in the LETKF.


    for iv , var in enumerate(obj.fields.keys()) :
        if iv == 0 :
           [nlev,nlat,nlon]=np.shape( obj.fields[var]['data'] )
           tmp_data =np.zeros(( nlev,nlat,nlon,nvar))
           tmp_az =np.zeros(( nlev,nlat,nlon,nvar))
           tmp_ra =np.zeros(( nlev,nlat,nlon,nvar))
           tmp_el =np.zeros(( nlev,nlat,nlon,nvar))
           tmp_n  =np.zeros(( nlev,nlat,nlon,nvar)).astype(int)

        tmp_data[:,:,:,iv] = obj.fields[var]['data']
        tmp_az  [:,:,:,iv] = obj.fields[var]['az']
        tmp_el  [:,:,:,iv] = obj.fields[var]['el']
        tmp_ra  [:,:,:,iv] = obj.fields[var]['ra']
        tmp_n   [:,:,:,iv] = obj.fields[var]['nobs']
        tmp_error     [iv] = obj.fields[var]['error']
        tmp_id        [iv] = obj.fields[var]['id']

    #Filter grid points in which the number of data points is less than min_n observations
    min_n = 10  #TODO this should became an input parameter
    tmp_n[ tmp_n < min_n ]=0

    cs.write_radar(nlon=nlon,nlat=nlat,nlev=nlev,nvar=nvar,
                   data_in=tmp_data,ndata_in=tmp_n,
                   grid_az=tmp_az,grid_el=tmp_el,grid_ra=tmp_ra,
                   error=tmp_error,ido=tmp_id,lambdar=tmp_lambda,  
                   filename=filename,
                   radar_lon=obj.radar.longitude['data'],
                   radar_lat=obj.radar.latitude['data'] ,
                   radar_z=obj.radar.altitude['data'] )

    #Write summary.
    for iv , var in enumerate(obj.fields.keys()) :
        totalnobs=np.sum( tmp_n[:,:,:,iv] > 0 )
        if totalnobs > 0 :
           print('Total observations for obstype ' + str(np.array(obj.fields[var]['id'])) + ' = ' + str(totalnobs) )
           print('Max value is ' + str(np.max( tmp_data[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ] )) )
           print('Min value is ' + str(np.min( tmp_data[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ] )) )
           print('Grid properties')
           print('Max azimuth is  ',np.max(tmp_az[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ])) 
           print('Min azimuth is  ',np.min(tmp_az[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ]))
           print('Max elevation is  ',np.max(tmp_el[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ]))
           print('Min elevation is  ',np.min(tmp_el[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ]))
           print('Max range is  ',np.max(tmp_ra[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ]))
           print('Min range is  ',np.min(tmp_ra[:,:,:,iv][ tmp_n[:,:,:,iv] > 0 ]))


def str2date(string):
    """ String must be with format %Y-%m-%dT%H:%M:%SZ """
    return datetime.strptime(string, '%Y-%m-%dT%H:%M:%SZ')

def date2str(date):
    """ Datetime object to string with format %Y-%m-%dT%H:%M:%SZ """
    return datetime.strftime(date, '%Y%m%d%H%M%S')

def datespan(startDate, endDate, delta):
    '''
    La funcion devuelve un "generator" que contiene un objecto date
    Input:
        starDate (objeto): de la clase datetime que indica la fecha inicial
        endDate (objeto): de la clase datetime que indica la fecha final
        delta (objeto): de la clase datetime que indica el intervalo temporal
    '''
    currentDate = startDate
    while currentDate < endDate:
        yield currentDate
        currentDate += delta

if __name__ == '__main__':
    file = 'cfrad.20091117_174348.000_to_20091117_174737.000_PAR_SUR.nc3'
    #file = './cfrad.20100111_000003.000_to_20100111_000340.001_ANG240_v1_SUR.nc'
    file = 'cfrad.20170926_171335.0000_to_20170926_171446.0000_RMA1_0122_03.nc3'
    #print(file)

    a = main_radar_so(file, 300, [2000, 2000, 25e3], {'ZH':[4001, 5, 0], 'VRAD':[4002, 2]})

    '''
    original = pyart.io.read(file)
    display = pyart.graph.RadarDisplay(original)
    print('PLOTTING ORIGINAL')
    for sweep in range(original.nsweeps):
        fig = plt.figure()
        display.plot_ppi('VRAD', sweep=sweep)
        plt.show()

    print('DOING SO')
    so = SoFields(file, ['VRAD'], [0, original.nrays],[10e3, 1000, 25e3])

    print('PLOTTING SO')
    for lev in range(so.grid.nlev):
        print(lev)
        fig = plt.figure()
        plt.pcolormesh(so.fields['grid_VRAD']['data'][lev,:,:])
        plt.colorbar()
        plt.show()
    '''



