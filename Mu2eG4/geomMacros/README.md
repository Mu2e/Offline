# Geometry Macros

This directory contains two useful ROOT macros that can be useful for inspecting the geometry

1. ```where_is_volume.C``` takes a string and will print all volumes that have that string in the name. It can take three arguments:
  * ```search_name``` is the name you are searching for (required)
  * ```gdmlname``` is the name of the GDML file (optional, default = ```mu2e.gdml```)
  * ```verbose``` is a boolean to print additional information about the volume (optional, default = ```false```)

Examples:
```
root -l -b -q where_is_volume.C\(\"VirtualDetector\"\)
root -l -b -q where_is_volume.C\(\"VirtualDetector\"\,\"mu2e.gdml\")
root -l -b -q where_is_volume.C\(\"VirtualDetector\"\,\"mu2e.gdml\",true)
```

1. 
