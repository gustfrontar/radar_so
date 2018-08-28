import numpy as np
import struct
import os

obsrecords = {
'speedy': {'names':['elm','lon','lat','lev','dat','err'],
           'formats':['i4','f4','f4','f4','f4','f4']},
'speedy2': {'names':['elm','lon','lat','lev','dat','err','ohx','oqc'],
           'formats':['i4','f4','f4','f4','f4','f4','f4','i4']},
'gfs':    {'names':['elm','lon','lat','lev','dat','err','typ'],
           'formats':['i4','f4','f4','f4','f4','f4','i4']},
'gfs2':   {'names':['elm','lon','lat','lev','dat','err','typ','tdf','ohx','oqc'],
           'formats':['i4','f4','f4','f4','f4','f4','i4','f4','f4','i4']},
'gfsdep': {'names':['elm','lon','lat','lev','dat','err','typ','tdf','omb','oma'],
           'formats':['i4','f4','f4','f4','f4','f4','i4','f4','f4','f4']},
'scale':  {'names':['elm','lon','lat','lev','dat','err','typ','tdf'],
           'formats':['i4','f4','f4','f4','f4','f4','i4','f4']},
'scale2': {'names':['set','idx','ohx','oqc','gri','grj'],
           'formats':['i4','i4','f4','i4','f4','f4']}
}

obselems = {
 2819: 'U',
 2820: 'V',
 3073: 'T',
 3074: 'Tv',
 3330: 'Q',
 3331: 'RH',
 4001: 'RADAR_REF',
 4002: 'RADAR_VR',
 4003: 'RADAR_PRH',
14593: 'PS',
19999: 'PRC',
99991: 'TCX',
99992: 'TCY',
99993: 'TCP'
}

obstypes = {
 0: 'OTHERS',
 1: 'ADPUPA',
 2: 'AIRCAR',
 3: 'AIRCFT',
 4: 'SATWND',
 5: 'PROFLR',
 6: 'VADWND',
 7: 'SATEMP',
 8: 'ADPSFC',
 9: 'SFCSHP',
10: 'SFCBOG',
11: 'SPSSMI',
12: 'SYNDAT',
13: 'ERS1DA',
14: 'GOESND',
15: 'QKSWND',
16: 'MSONET',
17: 'GPSIPW',
18: 'RASSDA',
19: 'WDSATR',
20: 'ASCATW',
21: 'TMPAPR',
22: 'PHARAD',
}

def endian_det(endian=None):
    if endian is None:
        return (np.dtype('f4'), 'i')
    elif endian == 'b' or endian == '>':
        return (np.dtype('>f4'), '>i')
    elif endian == 'l' or endian == '<':
        return (np.dtype('<f4'), '<i')
    else:
        raise ValueError("'endian' has to be [None|'b'|'l'|'>'|'<'].")

def readobs(fo, otype='gfs', endian=None):
    try:
        obsrecords[otype]
    except KeyError:
        raise ValueError("'otype' has to be " + str(list(obsrecords.keys())))
    reclen = len(obsrecords[otype]['names'])
    dtyp, dtypistr = endian_det(endian)

    if not fo.read(4):
        return False
    sglobs = np.fromfile(fo, dtype=dtyp, count=reclen).tolist()
    if not fo.read(4):
        return False

    for i in range(reclen):
        if obsrecords[otype]['formats'][i] == 'i4':
            sglobs[i] = int(round(sglobs[i]))
    return sglobs

def readobs_all(fo, otype='gfs', endian=None):
    try:
        obsrecords[otype]
    except KeyError:
        raise ValueError("'otype' has to be " + str(list(obsrecords.keys())))
    reclen = len(obsrecords[otype]['names'])
    dtyp, dtypistr = endian_det(endian)

    data = []
    fo.seek(0)
    while fo.read(4):
        sglobs = np.fromfile(fo, dtype=dtyp, count=reclen).tolist()
        for i in range(reclen):
            if obsrecords[otype]['formats'][i] == 'i4':
                sglobs[i] = int(round(sglobs[i]))
        data.append(tuple(sglobs))
        fo.read(4)

    return np.array(data, dtype=np.dtype(obsrecords[otype]))

def readobs_all_fast(fo, otype='gfs', endian=None):
    try:
        obsrecords[otype]
    except KeyError:
        raise ValueError("'otype' has to be " + str(list(obsrecords.keys())))
    reclen = len(obsrecords[otype]['names'])
    dtyp, dtypistr = endian_det(endian)

    filelen = os.fstat(fo.fileno()).st_size
    nobs = int(filelen / (4 * (reclen+2)))

    fo.seek(0)

    dataraw = (np.fromfile(fo, dtype=dtyp, count=nobs*(reclen+2)).reshape(nobs, reclen+2))[:,1:reclen+1]

    return np.array([tuple(i) for i in dataraw], dtype=np.dtype(obsrecords[otype]))

def readobs_radar_all(fo, otype='gfs', endian=None):
    try:
        obsrecords[otype]
    except KeyError:
        raise ValueError("'otype' has to be " + str(list(obsrecords.keys())))
    reclen = len(obsrecords[otype]['names'])
    dtyp, dtypistr = endian_det(endian)

    data = []
    fo.seek(0)

    fo.read(4)
    radar_lon = np.fromfile(fo, dtype=dtyp, count=1)[0]
    fo.read(4)
    fo.read(4)
    radar_lat = np.fromfile(fo, dtype=dtyp, count=1)[0]
    fo.read(4)
    fo.read(4)
    radar_z = np.fromfile(fo, dtype=dtyp, count=1)[0]
    fo.read(4)

    while fo.read(4):
        sglobs = np.fromfile(fo, dtype=dtyp, count=reclen).tolist()
        for i in range(reclen):
            if obsrecords[otype]['formats'][i] == 'i4':
                sglobs[i] = int(round(sglobs[i]))
        data.append(tuple(sglobs))
        fo.read(4)

    return radar_lon, radar_lat, radar_z, np.array(data, dtype=np.dtype(obsrecords[otype]))

def readobs_radar_all_fast(fo, otype='gfs', endian=None):
    try:
        obsrecords[otype]
    except KeyError:
        raise ValueError("'otype' has to be " + str(list(obsrecords.keys())))
    reclen = len(obsrecords[otype]['names'])
    dtyp, dtypistr = endian_det(endian)

    filelen = os.fstat(fo.fileno()).st_size
    nobs = int((filelen - 36) / (4 * (reclen+2)))

    fo.seek(0)

    fo.read(4)
    radar_lon = np.fromfile(fo, dtype=dtyp, count=1)[0]
    fo.read(4)
    fo.read(4)
    radar_lat = np.fromfile(fo, dtype=dtyp, count=1)[0]
    fo.read(4)
    fo.read(4)
    radar_z = np.fromfile(fo, dtype=dtyp, count=1)[0]
    fo.read(4)

    dataraw = (np.fromfile(fo, dtype=dtyp, count=nobs*(reclen+2)).reshape(nobs, reclen+2))[:,1:reclen+1]

    return radar_lon, radar_lat, radar_z, np.array([tuple(i) for i in dataraw], dtype=np.dtype(obsrecords[otype]))

def writeobs(fo, sglobs, endian=None):
    dtyp, dtypistr = endian_det(endian)

    sglobs_ndarray = np.array(tuple(sglobs), dtype=dtyp)
    fo.write(struct.pack(dtypistr, sglobs_ndarray.nbytes))
    sglobs_ndarray.tofile(fo)
    fo.write(struct.pack(dtypistr, sglobs_ndarray.nbytes))

def writeobs_all(fo, data, endian=None):
    dtyp, dtypistr = endian_det(endian)

    fo.seek(0)
    for sglobs in data:
        sglobs_ndarray = np.array(tuple(sglobs), dtype=dtyp)
        fo.write(struct.pack(dtypistr, sglobs_ndarray.nbytes))
        sglobs_ndarray.tofile(fo)
        fo.write(struct.pack(dtypistr, sglobs_ndarray.nbytes))

def writeobs_radar_all(fo, radar_lon, radar_lat, radar_z, data, endian=None):
    dtyp, dtypistr = endian_det(endian)

    fo.seek(0)

    tmp_ndarray = np.array(radar_lon, dtype=dtyp)
    fo.write(struct.pack(dtypistr, tmp_ndarray.nbytes))
    tmp_ndarray.tofile(fo)
    fo.write(struct.pack(dtypistr, tmp_ndarray.nbytes))
    tmp_ndarray = np.array(radar_lat, dtype=dtyp)
    fo.write(struct.pack(dtypistr, tmp_ndarray.nbytes))
    tmp_ndarray.tofile(fo)
    fo.write(struct.pack(dtypistr, tmp_ndarray.nbytes))
    tmp_ndarray = np.array(radar_z, dtype=dtyp)
    fo.write(struct.pack(dtypistr, tmp_ndarray.nbytes))
    tmp_ndarray.tofile(fo)
    fo.write(struct.pack(dtypistr, tmp_ndarray.nbytes))

    for sglobs in data:
        sglobs_ndarray = np.array(tuple(sglobs), dtype=dtyp)
        fo.write(struct.pack(dtypistr, sglobs_ndarray.nbytes))
        sglobs_ndarray.tofile(fo)
        fo.write(struct.pack(dtypistr, sglobs_ndarray.nbytes))
