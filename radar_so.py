"""
Functions to perform radar superobbing for LETKF-WRF data assimilation system
"""
from collections import defaultdict
import numpy as np
import numpy.ma
import pyart
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import manage_dates as dates

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
    def __init__(self, filename, variables, rays, grid):
        print('SO Fields object')
        self.radar = SoRadar(filename, variables, rays)
        self.grid = SoGrid(self.radar, grid[0], grid[1], grid[2])
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
            print(var)
            if self._check_field_exists(var):
                # Allocate vars in fields attribute
                self.allocate_var('grid_' + var, self.grid.nlon, self.grid.nlat, self.grid.nlev)
                # Get data from radar object
                vars_data[var] = self.radar.fields[var]

                # Convert dBZ to power
                if 'reflectivity' in vars_data[var]['standard_name']:
                    print('ToPower')
                    vars_data[var]['data'] = np.power(10, vars_data[var]['data']/10.)

        # Average data
        self.compute_grid_boxmean(vars_data)

    def allocate_var(self, var, nx, ny, nz):
        """ Allocate variables to a key in self.fields """
        self.fields[var]['nobs'] = np.zeros((nz, ny, nx))
        self.fields[var]['data'] = np.zeros((nz, ny, nx))
        self.fields[var]['az'] = np.zeros((nz, ny, nx))
        self.fields[var]['el'] = np.zeros((nz, ny, nx))
        self.fields[var]['ra'] = np.zeros((nz, ny, nx))

    def compute_grid_boxmean(self, variables):
        """
        Compute the i, j, k for each radar grid point and box-average the data.

        Parameters
        ----------
        variables : list of numpy array
        """
        print('BOX AVERAGE')
        for ia in range(self.radar.nrays):
            for ir in range(self.radar.ngates):

                # Get i, j, k using a very simple approach since we are assuming
                # a regular lat/lon/z grid
                [k, j, i] = self._radar2grid(ia, ir)

                # Skip data outside the grid domain
                if self._check_point_in_grid(i, j, k):
                    for var in variables.keys():
                        if not np.ma.is_masked(variables[var]['data'][ia, ir]):
                            #print(variables[var]['data'][ia, ir])
                            self.fields['grid_' + var]['data'][k, j, i] += \
                                variables[var]['data'][ia, ir]
                            #print(self.fields['grid_' + var]['data'][i, j, k])
                            self.fields['grid_' + var]['nobs'][k, j, i] += 1
                            self.fields['grid_' + var]['az'][k, j, i] += \
                                self.radar.azimuth['data'][ia]
                            self.fields['grid_' + var]['el'][k, j, i] += \
                                self.radar.elevation['data'][ia]
                            self.fields['grid_' + var]['ra'][k, j, i] += \
                                self.radar.range['data'][ir]

        # Compute observations boxaverage
        for var in variables.keys():
            nobs = self.fields['grid_' + var]['nobs']
            for key in self.fields['grid_' + var].keys():
                if not key == 'nobs':
                    self.fields['grid_' + var][key][nobs > 0] = \
                        self.fields['grid_' + var][key][nobs > 0]/nobs[nobs > 0]

            # Power to DBZ
            if var == 'TH' or var == 'dBZ':
                print('ToDBZ')
                tmp = self.fields['grid_' + var]['data']
                tmp[tmp > 0] = 10*np.log10(tmp[tmp > 0])
                tmp[tmp <= 0] = 0

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
        return True if 0 <= i < self.grid.nlon and  \
            0 <= j < self.grid.nlat and \
            0 <= k < self.grid.nlev else False

    #def _check_radarpoint_is_masked(self, point):
        #""" Check if a radar point has a nan value """
        #return True if np.ma.is_masked(point) else False

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
    radar : Radar
        PyArt-Like radar object
    """
    def __init__(self, filename, variables, ray_interval):

        self.__radar = pyart.io.read(filename)
        self.__variables = variables
        self.__ray_interval = ray_interval

        self.gate_longitude = {}
        self.gate_latitude = {}
        self.gate_altitude = {}
        self.azimuth = {}
        self.elevation = {}
        self.fields = defaultdict(dict)

        self.ngates = self.__radar.ngates
        self.longitude = self.__radar.longitude
        self.latitude = self.__radar.latitude
        self.nrays = ray_interval[-1] - ray_interval[0]

        self.get_data()

    def get_data(self):
        self.gate_longitude['data'] = self._extract_rays(self.__radar.gate_longitude['data'])
        self.gate_latitude['data'] = self._extract_rays(self.__radar.gate_latitude['data'])
        self.gate_altitude['data'] = self._extract_rays(self.__radar.gate_altitude['data'])
        self.azimuth['data'] = self._extract_rays(self.__radar.azimuth['data'])
        self.elevation['data'] = self._extract_rays(self.__radar.elevation['data'])
        self.range = self.__radar.range
        for var in self.__variables:
            self.fields[var] = self.__radar.fields[var]
            self.fields[var]['data'] = self._extract_rays(self.__radar.fields[var]['data'])

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
    def __init__(self, radar, dx, dz, maxz):

        self.dx = dx
        self.dz = dz

        # Compute possible value for `nlon` in order to cover the maximum radar range
        maxrange = 2.*np.max(radar.range['data'])
        self.nlon = np.ceil(maxrange/self.dx).astype('int')
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

def main_radar_so(filename, output_freq, grid_dims):
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
    """
    # Read radar volume using pyart
    radar = pyart.io.read(filename)

    # Get reflectivity and Doppler velocity variables names
    vars_name = []
    for key in radar.fields.keys():
        if key != 'CM' and key != 'TV' and ('reflectivity' in radar.fields[key]['standard_name'] or \
            'radial_velocity' in radar.fields[key]['standard_name']):
            vars_name.append(key)
    print(vars_name)

    # Get radar start and end times
    #TEMPORAL (ver si se puede obtener de algun atributo de radar en vez de usar
    #el nombre del archivo)
    #initime =  radar.metadata['time_coverage_start']
    #endtime =  radar.metadata['time_coverage_end']
    initime, endtime = parse_dates(filename)
    inidate = dates.str2date(initime)
    enddate = dates.str2date(endtime)

    # Create output files list according to out_freq
    output_dates = get_letkf_outputs(inidate, enddate, output_freq)

    #Loop over output files
    iray = 0
    inirayidx = 0
    for date in output_dates.keys():

        if check_file_exists(dates.date2str(date) + '.npz'):
            # Load obs and nobs
            pass

        # Get radar rays that contribute to current date
        diff = (output_dates[date][1]-inidate).total_seconds()
        while iray < radar.nrays and np.abs(radar.time['data'][iray] - diff) > 1e-3:
            iray += 1
        endrayidx = iray
        ray_limits = [inirayidx, endrayidx]
        print(ray_limits)

        # Initialize superobbing radar object
        #print(radar.fields[vars_name[0]]['data'].shape)
        #tmpradar = SoRadar(filename, vars_name, [inirayidx, endrayidx])
        #print(tmpradar.fields[vars_name[0]]['data'].shape)

        # Compute superobbing
        so = SoFields(filename, vars_name, ray_limits, grid_dims)

        # Write intermediate file

        # Write LETKF file

        inirayidx = iray

    return so


def check_file_exists(file):
    """ file : npz file """
    import os
    return True if os.path.exists(file) else False

def get_letkf_outputs(inidate, enddate, freq):
    """
    Create a list of letkf output filenames
    initime, endtime : datetime objects
    out_freq : int in seconds
    """
    if freq < 3600:
        # Get first possible output date
        inirounddate = datetime(inidate.year, inidate.month, inidate.day, inidate.hour, 0, 0)
        freqtimes = round((inidate-inirounddate).total_seconds()/freq)
        outputini = inirounddate + timedelta(seconds=freqtimes*freq)

        # Get last possible output date
        endrounddate = datetime(inidate.year, inidate.month, inidate.day, inidate.hour+1, 0, 0)
        freqtimes = round((endrounddate-enddate).total_seconds()/freq)
        outputend = endrounddate - timedelta(seconds=freqtimes*freq)

        # Create output list
        delta = timedelta(seconds=freq)
        output_dates = defaultdict(list)
        for date in dates.datespan(outputini, outputend+delta, delta):
            iniinterval = date - timedelta(seconds=freq/2)
            endinterval = date + timedelta(seconds=freq/2)
            output_dates[date] = [iniinterval, endinterval]

        # Check radar intial and final dates are in output list
        check_date_in_interval(inidate, output_dates[outputini][0], output_dates[outputini][1])
        check_date_in_interval(enddate, output_dates[outputend][0], output_dates[outputend][1])

        return output_dates
    else:
        raise ValueError('Not a valid output frequency:' + freq)

def check_date_in_interval(date, lower, upper):
    if not lower <= date <= upper:
        raise ValueError('Date not in interval: ' + date)
    return True

def write_grid(radar, so_ref, so_dv, grid, file):
    pass

def write_letkf(so_ref, so_dv, grid, file):
    pass

def parse_dates(filename):
    """ Get initial and final times from a radar filename """
    ini = filename.split('_')[0].split('.')[1] + filename.split('_')[1].split('.')[0]
    end = filename.split('_')[0].split('.')[1] + filename.split('_')[4].split('.')[0]

    #Si llegara a ser necesario discriminar entre radares
    #if 'PAR' in filename:
        # cfrad.20091117_174348.000_to_20091117_174737.000_PAR_SUR.nc3
    #elif 'ANG' in filename:
        # cfrad.20100111_000003.000_to_20100111_000340.001_ANG240_v1_SUR.nc
    #elif 'PER' in filename:
    #elif 'RMA1' in filename:
        # cfrad.20170926_171335.0000_to_20170926_171446.0000_RMA1_0122_03.nc3
    #else:
    #    raise ValueError('Radar not recognized')

    return ini, end

if __name__ == '__main__':
    file = 'cfrad.20091117_174348.000_to_20091117_174737.000_PAR_SUR.nc3'
    #file = './cfrad.20100111_000003.000_to_20100111_000340.001_ANG240_v1_SUR.nc'
    file = 'cfrad.20170926_171335.0000_to_20170926_171446.0000_RMA1_0122_03.nc3'
    print(file)

    #a = main_radar_so(file, 300, [2000, 2000, 25e3])

    #print('PLOTTING')
    original = pyart.io.read(file)
    display = pyart.graph.RadarDisplay(original)
    print('PLOTTING ORIGINAL')
    for i in range(original.nsweeps):
        fig = plt.figure()
        display.plot_ppi('VRAD', sweep=i)
        plt.show()

    print('DOING SO')
    so = SoFields(file, ['VRAD'], [0, original.nrays],[10e3, 1000, 25e3])

    print('PLOTTING SO')
    for i in range(so.grid.nlev):
        print(i)
        fig = plt.figure()
        plt.pcolormesh(so.fields['grid_VRAD']['data'][i,:,:])
        plt.colorbar()
        plt.show()
