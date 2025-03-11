#!/usr/bin/env python3
import numpy as np
from array import array

def read_ascii(filename):
    # read ASCII SWMF file "filename" and return a dictionary with the content
    f = open(filename,'r')
    headline = f.readline()
    step, time, ndim, nparam, nvar = f.readline().split()
    step = np.int32(step)
    time = np.float32(time)
    ndim = np.int32(ndim)
    nparam = np.int32(nparam)
    nvar = np.int32(nvar)
    size = np.array(f.readline().split(), dtype=np.int32)
    if nparam > 0:
        param = np.array(f.readline().split(), dtype=np.float32)
    namevar = f.readline()
    ngrid = np.prod(size)
    # read the remaining numbers (this will not work for .outs files)
    griddata = np.array(f.read().split(), dtype=np.float32)
    # reshape and transpose array to match default Python ordering
    griddata = griddata.reshape(ngrid, ndim+nvar).T
    f.close()
    # extract coordinates and variables
    coord = griddata[:ndim,:]
    var = griddata[ndim:,:]

    # debug
    #print('header=',headline)
    #print('step,time,ndim,nparam,nvar=', step, time, ndim, nparam,nvar)
    #print('size=',size)
    #print('ngrid=',ngrid)
    #if nparam > 0:
    #    print('param=', param)
    #print('namevar=', namevar)
    #print('coord=',coord)
    #print('var=',var)

    # Create dictionary
    data ={"headline": headline,
           "step"   : step,
           "time"   : time,
           "ndim"   : ndim,
           "size"   : size,
           "nparam" : nparam,
           "nvar"   : nvar,
           "namevar": namevar,
           "coord"  : coord,
           "var"    : var}
    if nparam > 0:
        data["param"] = param 
    return data
###############################################################################
def fortran_string(string, l):
    # create bytestring of length l with spaces added at the end
    if len(string) > l:
        print("ERROR in Fortran string:")
        print("Length of string=",len(string),"> l=",l)
        exit(1)
    return string.encode()[:-1] + b' '*(l+1 - len(string))
###############################################################################
def fortran_record(bytearray):
    # write a Fortran record with 4-byte markers at both ends
    reclen = np.int32(len(bytearray)).tobytes()
    return reclen + bytearray + reclen
###############################################################################
def write_real4(filename, data):
    # Print data in real4 format into file filename
    # The data dictionary must contain "size", "namevar", "coord", "var"
    # containing the grid size, coordinate+variable+param names,
    # coordinate array and variable array.
    # The data["param"] containing scalar parameters is optional.
    # All arrays should be 4-byte numpy arrays.

    f = open(filename,'wb')

    headline = data["headline"] if "headline" in data else "Missing headline"
    size = np.int32(data["size"])
    namevar = data["namevar"]
    step = np.int32(data["step"] if "step" in data else 0)
    time = np.float32(data["time"] if "time" in data else 0.0)
    param = data["param"] if "param" in data else []
    coord = data["coord"]
    var = data["var"]

    ndim = np.int32(len(size))
    nparam = np.int32(len(param))
    nvar = np.int32(len(namevar.split()) - ndim - nparam)

    # Check consistency
    ngrid = np.prod(size)
    if coord.size != ngrid*ndim:
        print("ERROR in write_real4: incorrect size of coord=", coord.size)
        print('size=', size, 'ndim=', ndim,' ngrid=',ngrid)
        exit(1)

    if var.size != ngrid*nvar:
        print("ERROR in write_real4: incorrect size of var=", var.size)
        print('size=', size, 'ndim=', ndim,' ngrid=', ngrid)
        print("namevar=", namevar," nparam=", nparam)
        exit(1)
        
    #print("write: nvar=", nvar)
    
    stringlength = 79 if len(headline) < 80 and len(namevar) < 80 else 500
    f.write(fortran_record(fortran_string(headline, stringlength)))
    f.write(fortran_record(step.tobytes()
                           + time.tobytes()
                           + ndim.tobytes()
                           + nparam.tobytes()
                           + nvar.tobytes()))
    f.write(fortran_record(size.tobytes()))
    if nparam > 0:
        f.write(fortran_record(param.tobytes()))
    f.write(fortran_record(fortran_string(namevar, stringlength)))
    f.write(fortran_record(coord.tobytes()))
    for i in range(nvar):
        f.write(fortran_record(var[i,:].tobytes()))
    f.close()
###############################################################################
# Main program with hard coded filenames for testing
filein  = 'dbH_plotfile.out'
fileout = 'dBH_real4.out'
print('reading   ', filein)
data = read_ascii(filein)
write_real4(fileout, data)
print('wrote out ', fileout)
