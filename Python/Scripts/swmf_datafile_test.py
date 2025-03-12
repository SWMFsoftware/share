#!/usr/bin/env python3
import swmf_datafile as swmf

# minimal data
data = {
    "name": "longitude latitude dBN dBE",
    "coord": [[[1.,2.,3.,4.],[1.,2.,3.,4.],[1.,2.,3.,4.]],
              [[1.,1.,1.,1.],[2.,2.,2.,2.],[3.,3.,3.,3.]]],
    "state": [[[10.,11.,12.,13.],[100.,110.,120.,130.],[1000.,1100.,1200.,1300]],
              [[20.,21.,22.,23.],[200.,210.,220.,230.],[2000.,2100.,2200.,2300]]]
}
swmf.show_data(data)

# write ASCII file with 3 decimals
fileout = 'file_ascii.outs'
swmf.write_file(data, fileout, format='10.3f')
# Add another snapshot
data["time"] = 1.0
swmf.write_file(data, fileout, format='10.3f', append=True)
print('wrote 2 snapshots into ', fileout)
print("Type of "+fileout+":", swmf.file_format(fileout))

# write single precision binary file
fileout = 'file_real4.out'
swmf.write_file(data, fileout, "real4")
print("wrote out ", fileout)
print("Type of "+fileout+":", swmf.file_format(fileout))

# Add head line and dimensions
data["head"] = "Proper header"
data["cart"] = False
data["dims"] = [1, 12] # typical for AMR grid
data["pars"] = [0.1, 0.2, 0.3, 0.4, 0.5] # add 5 scalar parameters
# New data name is more than 79 characters
data["name"] += " ParameterName1 ParameterName2 ParameterName3 ParameterName4 ParameterName5"

# write double precision binary file
fileout = 'file_real8.outs'
swmf.write_file(data, fileout, "real8")
# add another snapshot
data["time"] = 2.0
swmf.write_file(data, fileout, "real8", append=True)
print('wrote 2 snapshots into ', fileout)
print("Type of "+fileout+":", swmf.file_format(fileout))

# read back ascii file
filein = 'file_ascii.outs'
print('reading  ', filein)
data2 = swmf.read_file(filein, verbose=True)
print("head=",data2["head"])
print("dims=",data2["dims"], "cart=", data2["cart"])
print("name=",data2["name"])

# read back single precision binary file
filein = 'file_real4.out'
print('reading  ', filein)
data2 = swmf.read_file(filein)
swmf.show_data(data2)
    
# read back double precision binary file
filein = 'file_real8.outs'
f = open(filein,'rb')
print('reading  ', filein)
data2 = swmf.read_file(f, verbose=True)
print('read next', filein)
data3 = swmf.read_file(f, verbose=True)
f.close
