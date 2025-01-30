//
// Purpose:
// --------
// This function identifies and visualizes volumes intersecting a 3D bounding
// box in a ROOT-based geometry using GDML files.
// It maps user-specified coordinates (based on reference drawings) into the
// internal coordinate system and performs a grid-based search within the
// bounding box.
//
// Inputs:
// -------
//  - Bounding Box Center Coordinates:
//      - X Coordinate (x0):
//        Interpolated floating-point value between reference lines in the
//        drawing.
//        Unit lengths are:
//              40' for locations north of reference line 2 (values < 2).
//              25' for locations south of reference line 2 (values > 2).
//      - Y Coordinate (y0):
//        Elevation quoted in feet based on drawings.
//        Use fractions of feet, avoiding inches.
//        Reference elevation: 746.5' (746' 6").
//      - Z Coordinate (z0_refline):
//        Interpolated floating-point value between vertical lines
//        (A, B, C... → 1, 2, 3...).
//        Unit lengths are:
//              27' for locations east of reference line B (values < 2).
//              26' for locations west of reference line B (values > 2).
//
// - Bounding Box Dimensions
//      Half-Lengths of Box Sides (halfdx_feet, halfdy_feet, halfdz_feet):
//      Specify the dimensions of the box in feet;
//      converted internally to millimeters.
//
// - GDML File Prefix
//      The prefix of the GDML file containing the geometry.
//      Example: "mu2e_40".
//
// - Step Fraction
//      Fraction of box dimensions used as the grid step size during the search.
//      Smaller values provide finer resolution but may increase runtime.
//
// Key Features:
// -------------
//
// - Coordinate System Mapping:
//      Translates user-specified reference coordinates (x0, y0, z0) into the
//      internal geometry coordinate system (millimeters).
//      Uses interpolation rules derived from the provided reference drawings.
//      (doc-db=11669 filename:SC-SET_6-10-2-AB.pdf)
//
// - Bounding Box Dimensions
//      Defines a 3D bounding box centered on the computed coordinates, with
//      dimensions based on the provided half-lengths in feet (converted to mm)
//
// - Volume Search
//      A grid-based search is performed within the box, with step sizes
//      determined by stepFraction.
//      At each grid point, the macro identifies the intersecting volume using
//      the find_volume_at_point function.
//      Logs volume names and their intersection coordinates.
//
// - Volume Visualization
//      Initially hides all volumes for clarity.
//      Makes intersecting volumes visible and assigns distinct colors to the
//      closest ones for easy differentiation.
//      Visualized in the ROOT OpenGL viewer for interactive exploration.
//
// Output:
// -------
//  - Logs the names and coordinates of intersecting volumes.
//  - Opens a TBrowser for interactive geometry exploration.
//  - Allows clipping plane manipulation to inspect volume boundaries.
//
// How to Run:
// -----------
//
// - Determine the box center coordinates using the reference drawings.
//   Example:
//       X Coordinate: A bit north of reference line 2 → 1.8.
//       Y Coordinate: Elevation 743' → 743.0.
//       Z Coordinate: A bit east of reference line A (corresponds to 1) → 0.9.
//
// - Specify bounding box dimensions (in feet) and the GDML file prefix.
//
// - Execute the macro with the desired parameters.
//   Example:
//
//   root -l explore_building_area.C\(1.8,743.,0.9,2,2,2,\"mu2e_40\",0.01\)
//
//        1.8, 743., 0.9: Box center coordinates.
//        2, 2, 2: Half-dimensions of the bounding box (in feet).
//        "mu2e_40": GDML file prefix.
//        0.01: Step fraction for grid resolution.
// Notes:
// ------
//  - Adjust the bounding box size and step fraction for optimal performance
//    and clarity.
//  - Use clipping planes in the ROOT viewer to inspect volume intersections
//    effectively.
//  - (R. Kutschke's tip) In case the OpenGL viewer is not working try in a
//    clean shell ( not even mu2einit ):
//  > source /cvmfs/scisoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
//  > spack load root/xarsy2v
//
// Author: S. Di Falco
// Jan 23, 2025

#include <iostream>
#include <fstream>
#include <vector>
#include "clean_volume_names.C"
#include "TGeoVolume.h"
#include "TSystem.h"

// Function to find the volume name at a given position
TString find_volume_at_point(double* pos, TGeoManager* geom) {
  double x = pos[0] / 10.; // Convert mm to cm
  double y = pos[1] / 10.;
  double z = pos[2] / 10.;
  TGeoNode* node = geom->FindNode(x, y, z);

  if (!node) {
    std::cerr << "The position (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ") is outside the world.\n";
    return "";
  }
  return clean_volume_names(node->GetName());
}

// Function to visualize volumes within a specified bounding box
int explore_building_area(
  double x0 = 0.,
  double y0 = 0.,
  double z0 = 0.,
  double halfdx_feet = 5.,
  double halfdy_feet = 5.,
  double halfdz_feet = 5.,
  TString gdmlFilePrefix = "mu2e_40",
  double stepFraction = 0.05) {

  // Conversion factor: feet to mm
  const double feet2mm = 304.8;

  // ********** Load Geometry **********
  TString gdmlFile = gdmlFilePrefix + ".gdml";
  if (gSystem->AccessPathName(gdmlFile.Data())) {
    std::cerr << "Error: GDML file not found: " << gdmlFile << std::endl;
    return -1;
  }

  TGeoManager* geom = TGeoManager::Import(gdmlFile);
  if (!geom) {
    std::cerr << "Error: Failed to load geometry from " << gdmlFile << std::endl;
    return -2;
  }

  // ********** Define Bounding Box **********
  double bbox_origin[3];

  // X Coordinate Conversion
  const double x_base = -9212.6;     // mm
  const double xstep_north = 12192.; // mm (40')
  const double xstep_south = 7620.;  // mm (25')
  bbox_origin[0] = (x0 <= 2.)
    ? x_base + xstep_north * (2. - x0)
    : x_base - xstep_south * (x0 - 2.);

  // Y Coordinate Conversion
  const double y_main_mm = 5460.4;  // mm
  const double y_main_feet = 746.5; // feet
  bbox_origin[1] = y_main_mm + (y0 - y_main_feet) * feet2mm;

  // Z Coordinate Conversion
  const double z_base = 2133.6;       // mm
  const double zstep_east = 8229.6;   // mm (27')
  const double zstep_west = 7924.8;   // mm (26')
  bbox_origin[2] = (z0 <= 2.)
    ? z_base - zstep_east * (2. - z0)
    : z_base + zstep_west * (z0 - 2.);

  // Box Dimensions
  double halfLength[3] = {
    halfdx_feet * feet2mm,
    halfdy_feet * feet2mm,
    halfdz_feet * feet2mm
  };

  // Configure Steps
  int nsteps = int(2. / stepFraction);
  double step[3];
  for (int i = 0; i < 3; i++) {
    step[i] = 2. * halfLength[i] / nsteps;
  }

  // ********** Volume Search **********
  std::vector<TString> volnames;

  // Hide all volumes initially
  for (const auto obj : *geom->GetListOfVolumes()) {
    TGeoVolume* vol = (TGeoVolume*)obj;
    vol->SetVisibility(0);
  }

  // Define colors for visualization
  std::vector<int> colors = {kOrange, kRed, kBlue, kMagenta, kCyan, kSpring};
  int coloredVolumeCount = 0;

  // Grid Scan to find intersecting volumes
  double gridPoint[3]; // Coordinates in mm
  for (int ix = 0; ix < nsteps; ++ix) {
    gridPoint[0] = bbox_origin[0] + ((ix % 2 == 0) ? ix / 2 * step[0] : -(ix + 1) / 2 * step[0]);;

    for (int iy = 0; iy < nsteps; ++iy) {
      gridPoint[1] = bbox_origin[1] + ((iy % 2 == 0) ? iy / 2 * step[1] : -(iy + 1) / 2 * step[1]);

      for (int iz = 0; iz < nsteps; ++iz) {
        gridPoint[2] = bbox_origin[2] + ((iz % 2 == 0) ? iz / 2 * step[2] : -(iz + 1) / 2 * step[2]);

        TString vol_name = find_volume_at_point(gridPoint, geom);
        if (vol_name != "HallAir" && !vol_name.IsNull()) {
          if (std::find(volnames.begin(), volnames.end(), vol_name) == volnames.end()) {
            TGeoVolume* vol = geom->FindVolumeFast(vol_name.Data());
            if (vol) {
              volnames.push_back(vol_name);
              vol->SetVisibility(1);
              if (coloredVolumeCount < colors.size()) {
                vol->SetLineColor(colors[coloredVolumeCount]);
                ++coloredVolumeCount;
              }
            }
            std::cout << "Volume: " << vol_name << " at (" << gridPoint[0]
                      << ", " << gridPoint[1] << ", " << gridPoint[2] << ") mm\n";
          }
        }
      }
    }
  }

  // Final output and visualization
  std::cout << volnames.size() << " volumes intersecting the bounding box.\n";
  TBrowser* browser = new TBrowser();
  geom->GetVolume("HallAir")->Draw("ogl");

  return 0;
}
