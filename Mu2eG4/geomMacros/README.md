# Geometry Macros

This directory contains two useful ROOT macros that can be useful for inspecting the geometry

## ```where_is_volume.C```
This macro takes a string and will print all volumes that have that string in the name. It can take three arguments:
* ```search_name``` is the name you are searching for (required)
* ```gdmlname``` is the name of the GDML file you want to inspect (optional, default = ```mu2e.gdml```)
* ```verbose``` is a boolean to print additional information about the volume (optional, default = ```false```)

Examples:
```
root -l -b -q where_is_volume.C\(\"VirtualDetector\"\)
root -l -b -q where_is_volume.C\(\"VirtualDetector\"\,\"mu2e.gdml\")
root -l -b -q where_is_volume.C\(\"VirtualDetector\"\,\"mu2e.gdml\",true)
```

## ```find_volume_at_point.C```
This macro takes an x, y, z position (in Mu2e global coordinates) and prints information for the volume at that position. Some warnings:
* this will return the deepest volume in the hierarchy, and
* if two volumes are touching at that point, the volume returned is ambiguous

There are two functions defined in this macro. One can be used on the command line and takes the following arguments:
* ```x, y, z``` is the position (required)
* ```gdmlname``` is the name of the GDML file you want to inspect (optional, default = ```mu2e.gdml```)

Examples:
```
root -l -b -q find_volume_at_point.C\(-3904,0,4000\)
root -l -b -q find_volume_at_point.C\(-3904,0,4000,\"mu2e.gdml\"\)
````

The other can be called from another macro so that you can ask for multiple points without having to reload the GDML each time. This takes the following arguments:
* ```x, y, z``` is the position (required)
* ```geom``` is the ```TGeoManager``` (required)

Example: see ```find_stm_beamline_volumes.C```
