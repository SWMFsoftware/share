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

print("name=", data["name"])
print("coord.shape=", data["coord"])
print("state.shape=", data["state"])

# write ASCII file with 3 decimals
fileout = 'file_ascii.out'
swmf.write_file(data, fileout, '10.3f')
print("Type of "+fileout+":", swmf.file_format(fileout))
    
# write single precision binary file
fileout = 'file_real4.out'
swmf.write_file(data, fileout, "real8")
print('wrote out ', fileout)
print("Type of "+fileout+":", swmf.file_format(fileout))

# Add head line and dimensions
data["head"] = "Proper header"
data["cart"] = False
data["dims"] = [1, 12] # typical for AMR grid
data["pars"] = [0.1, 0.2, 0.3, 0.4, 0.5] # add 5 scalar parameters
# New data name is more than 79 characters
data["name"] += " ParameterName1 ParameterName2 ParameterName3 ParameterName4 ParameterName5"

# write double precision binary file
fileout = 'file_real8.out'
swmf.write_file(data, fileout, "real8")
print('wrote out ', fileout)
print("Type of "+fileout+":", swmf.file_format(fileout))

# read back ascii file
filein = 'file_ascii.out'
print('reading  ', filein)
data2 = swmf.read_file(filein)
print("head=",data2["head"])
print("dims=",data2["dims"], "cart=", data2["cart"])
print("name=",data2["name"])

# read back single precision binary file
filein = 'file_real4.out'
print('reading  ', filein)
data2 = swmf.read_file(filein)
print("head=",data2["head"])
print("dims=",data2["dims"], "cart=", data2["cart"])
print("name=",data2["name"])
    
# read back double precision binary file
filein = 'file_real8.out'
print('reading  ', filein)
data2 = swmf.read_file(filein)
print("head=",data2["head"])
print("dims=",data2["dims"], "cart=", data2["cart"])
print("name=",data2["name"])
