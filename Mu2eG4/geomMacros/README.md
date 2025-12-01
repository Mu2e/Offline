# Geometry Macros

This directory contains some scripts that are useful for inspecting the geometry.

## ```where_is_volume.C```
This macro takes a string and will print all volumes that have that string in the name. It can take three arguments:
* ```search_name``` is the name you are searching for (required)
* ```gdmlname``` is the name of the GDML file you want to inspect (optional, default = ```mu2e_common.gdml```)
* ```verbose``` is a boolean to print additional information about the volume (optional, default = ```false```)

Examples:
```
root -l -b -q 'where_is_volume.C("VirtualDetector")'
root -l -b -q 'where_is_volume.C("VirtualDetector","mu2e_common.gdml")'
root -l -b -q 'where_is_volume.C("VirtualDetector","mu2e_common.gdml",true)'
```

## ```find_volume_at_point.C```
This macro takes an x, y, z position (in Mu2e global coordinates) and prints information for the volume at that position. Some warnings:
* this will return the deepest volume in the hierarchy
* if two volumes are touching at the given position, the volume returned is seemingly random
   * e.g. if volume A is fully contained within volume B, and you ask for a point on the boundary of volume A, round-off errors could give you volume B instead of the expected volume A
* ```StepPointMCs``` have positions on volume boundary

There are two functions defined in this macro. One can be used on the command line and takes the following arguments:
* ```x, y, z``` is the position (required)
* ```gdmlname``` is the name of the GDML file you want to inspect (optional, default = ```mu2e_common.gdml```)
* ```output_csv``` is a boolean for if you want the output printined in csv format (optional, default = ```false```)
   * note that this removes the memory addresses from the names of the volumes, which are used to distinguish different volumes that have the same name
* ```out``` is an ```std::ostream``` if you want to redirect the output to a file (optional, default = ```std::cout```)

Examples:
```
root -l -b -q 'find_volume_at_point.C(-3904,0,4000)'
root -l -b -q 'find_volume_at_point.C(-3904,0,4000,"mu2e_common.gdml")'
````

The other can be called from another macro so that you can ask for multiple points without having to reload the GDML each time. This takes the following arguments:
* ```x, y, z``` is the position (required)
* ```geom``` is the ```TGeoManager``` (required)
* ```output_csv``` is a boolean for if you want the output printined in csv format (optional, default = ```false```)
* ```out``` is an ```std::ostream``` if you want to redirect the output to a file (optional, default = ```std::cout```)x

Example: see ```find_stm_beamline_volumes.C``` and ```write_stm_beamline_volumes_to_csv.C```

## ```write_geom_diff_to_csv.py```
This python script takes two csv files created by ```find_volume_at_point```, merges them, and then writes the rows where there are differences in the volume name, material, or path to an new csv file.

Note: this requires you are in a python environment with the ```pandas``` package
