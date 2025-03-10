#!/usr/bin/env python3
import numpy as np

# Read and write SWMF plotfiles in ascii, real4 and real8 formats
# After reading the data dictionary contains
#   "head": headline string
#   "step": time step (integer)
#   "time": time (real)
#   "ndim": number of coordinates (1, 2, or 3)
#   "npar": number of scalar parameters (0, 1, ...)
#   "nvar": number of variables in the state (1, 2, ...)
#   "dims": array of ndim integers giving the grid dimension
#   "pars": array of npar reals with scalar parameters (if npar > 0)
#   "name": coordinate, variable and parameter names separated by spaces
#   "coord": coordinate array (ndim,ngrid)
#   "state": state array (nvar,ngrid)
#
# For writing the data dictionary must contain "name", "coord", "state",
# while "pars" is optional. If "dims" is present it is used and checked
# against the shape of "coord".

###############################################################################
def file_format(filename):
    # Return the file format
    f = open(filename, 'rb')
    # check if first four characters look like a record length
    reclen = int.from_bytes(f.read(4),'little')
    if reclen != 79 and reclen != 500:
        f.close()
        return "ascii"
    # skip the head line
    f.read(reclen+4)
    reclen = int.from_bytes(f.read(4),'little')
    f.close()
    if reclen == 20:
        return "real4"
    elif reclen == 24:
        return "real8"
    else:
        return "ascii"
    print('bytes=', bytes, 'reclen=', reclen)
    return "ascii"
    
###############################################################################
def read_plotfile(filename, fileformat='unknown'):
    if fileformat == "unknown":
        # figure out the file format
        fileformat = file_format(filename)

    if fileformat == "real4" or fileformat == "real8":
        return read_binary(filename)
    else:
        return read_ascii(filename)

###############################################################################
def write_plotfile(data, filename="plotfile.out", fileformat="18.10e"):
    if fileformat == 'real4' or fileformat == 'real8':
        write_binary(data, filename, fileformat)
    else:
        write_ascii(data, filename, fileformat)
    return

###############################################################################
def read_ascii(filename):
    # read ASCII SWMF file "filename" and return a dictionary with the content
    f = open(filename,'r')
    head = f.readline()
    step, time, ndim, npar, nvar = f.readline().split()
    step = np.int32(step)
    time = np.float64(time)
    ndim = np.int32(ndim)
    npar = np.int32(npar)
    nvar = np.int32(nvar)
    dims = np.array(f.readline().split(), dtype=np.int32)
    if npar > 0:
        pars = np.array(f.readline().split(), dtype=np.float64)
    name = f.readline()
    ngrid = np.prod(dims)
    # read the remaining numbers (this will not work for .outs files)
    griddata = np.array(f.read().split(), dtype=np.float64)
    # reshape and transpose array to match default Python ordering
    griddata = griddata.reshape(ngrid, ndim+nvar).T
    f.close()
    # extract coordinates and variables
    coord = griddata[:ndim,:]
    state = griddata[ndim:,:]

    # debug
    #print('head=',head)
    #print('step,time,ndim,npar,nvar=', step, time, ndim, npar, nvar)
    #print('dims=',dims)
    #print('ngrid=',ngrid)
    #if npar > 0:
    #    print('pars=', pars)
    #print('name=', name)
    #print('coord=',coord)
    #print('state=',state)

    # Create dictionary
    data ={"head"   : head,
           "step"   : step,
           "time"   : time,
           "ndim"   : ndim,
           "dims"   : dims,
           "npar"   : npar,
           "nvar"   : nvar,
           "name"   : name,
           "coord"  : coord,
           "state"  : state}
    if npar > 0:
        data["pars"] = pars 
    return data
###############################################################################
def read_binary(filename):
    return
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
def write_ascii(data, filename="plotfile.out", fileformat="18.10e"):
    # Print data into ascii file filename using fileformat for real arrays.

    if not "name" in data or not "state" in data or not "coord" in data:
        print("ERROR in write_ascii: missing name, state or coord in data")
        print("Could not write file", filename)
        return 1

    f = open(filename,'w')

    name = data["name"]
    coord = np.float64(data["coord"])
    state = np.float64(data["state"])
    head = data["head"] if "head" in data else "Missing head"
    dims = data["dims"] if "dims" in data else coord.shape[1:]
    step = data["step"] if "step" in data else 0
    time = data["time"] if "time" in data else 0.0
    pars = np.float64(data["pars"]) if "pars" in data else []

    ndim = len(dims)
    npar = len(pars)
    nvar = len(name.split()) - ndim - npar

    # Check consistency
    ngrid = np.prod(dims)
    if coord.size != ngrid*ndim:
        print("ERROR in write_ascii: incorrect size of coord=", coord.size)
        print('dims=', dims, 'ndim=', ndim,' ngrid=',ngrid)
        exit(1)

    if state.size != ngrid*nvar:
        print("ERROR in write_ascii: incorrect size of state=", state.size)
        print('dims=', dims, 'ndim=', ndim,' ngrid=', ngrid)
        print("name=", name," npar=", npar)
        exit(1)

    float_formatter = "{:"+fileformat+"}"

    print("float_formatter =",float_formatter)
    np.set_printoptions(formatter={'float_kind':float_formatter.format})
    f.writelines(head)
    template = "{0:d} {1:"+fileformat+"} {2:d} {3:d} {4:d}\n"
    f.write(template.format(step, time, ndim, npar, nvar))
    f.write(str(dims)[1:-1]+"\n")
    template = "{:18.10e}"
    if npar > 0:
        f.write(str(pars)[1:-1]+"\n")
    f.writelines(name)
    print("coord.shape=", coord.shape)
    for i in range(ngrid):
        f.write(str(coord[:,i])[1:-1] + str(state[:,i])[1:-1] + "\n")

    f.close()
###############################################################################
def write_binary(data, filename="plotfile.out", fileformat="real4"):
    # Print binary data in fileformat (real4 or real8) into file filename.
    # All arrays should be numpy arrays.

    print("starting write_binary: filename, fileformat=", filename, fileformat)
    
    if not "name" in data or not "state" in data or not "coord" in data:
        print("ERROR in write_plotfile: missing name, state or coord in data")
        print("Could not write file", filename)
        return 1

    f = open(filename,'wb')

    # required
    name = data["name"]
    coord = data["coord"]
    state = data["state"]
    time = data["time"] if "time" in data else 0.0
    pars = data["pars"] if "pars" in data else []
    npar = np.int32(len(pars))

    # convert to real4 / real8 as needed
    if fileformat == "real8":
        if coord.dtype != "float64":
            coord = np.float64(coord)
        if state.dtype != "float64":
            state = np.float64(state)
        time = np.float64(time)
        if npar > 0 and pars.dtype != "float64":
            pars = np.float64(pars)
            
    else:
        if coord.dtype != "float32":
            coord = np.float32(coord)
        if state.dtype != "float32":
            state = np.float32(state)
        time = np.float32(time)
        if npar > 0 and pars.dtype != "float32":
            pars = np.float32(pars)
            
    head = data["head"] if "head" in data else "Missing head"
    dims = np.int32(data["dims"]) if "dims" in data else coord.shape[1:]
    step = np.int32(data["step"] if "step" in data else 0)

    ndim = np.int32(len(dims))
    nvar = np.int32(len(name.split()) - ndim - npar)

    # Check consistency
    ngrid = np.prod(dims)
    if coord.size != ngrid*ndim:
        print("ERROR in write_binary: incorrect size of coord=", coord.size)
        print('dims=', dims, 'ndim=', ndim,' ngrid=',ngrid)
        exit(1)

    if state.size != ngrid*nvar:
        print("ERROR in write_binary: incorrect size of state=", state.size)
        print('dims=', dims, 'ndim=', ndim,' ngrid=', ngrid)
        print("name=", name," npar=", npar)
        exit(1)
        
    stringlength = 79 if len(head) < 80 and len(name) < 80 else 500
    f.write(fortran_record(fortran_string(head, stringlength)))
    f.write(fortran_record(step.tobytes()
                           + time.tobytes()
                           + ndim.tobytes()
                           + npar.tobytes()
                           + nvar.tobytes()))
    f.write(fortran_record(dims.tobytes()))
    if npar > 0:
        f.write(fortran_record(pars.tobytes()))
    f.write(fortran_record(fortran_string(name, stringlength)))
    f.write(fortran_record(coord.tobytes()))
    for i in range(nvar):
        f.write(fortran_record(state[i,:].tobytes()))
    f.close()
    return

###############################################################################
def test_plotfile():
    filein  = 'file_ascii.ini'
    fileout = 'file_real4.out'
    print("Type of "+filein+":", file_format(filein))
    print('reading   ', filein)
    data = read_plotfile(filein)
    write_plotfile(data, fileout, "real4")
    print('wrote out ', fileout)
    print("Type of "+fileout+":", file_format(fileout))

    fileout = 'file_real8.out'
    write_plotfile(data, fileout, "real8")
    print("Type of "+fileout+":", file_format(fileout))

    fileout = 'file_ascii.out'
    write_plotfile(data, fileout, '10.3f')
    print("Type of "+fileout+":", file_format(fileout))

    return
###############################################################################
test_plotfile()
