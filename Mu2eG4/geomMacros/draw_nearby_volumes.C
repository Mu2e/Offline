//
// Purpose:
// --------
// This ROOT Macro identifies and visualizes volumes in a ROOT geometry that
// are close to a specified volume.
//
// Inputs:
// -------
//
// - Volume Name (TString):
//     The name of the volume to inspect and highlight.
//
// - Margin (double, in mm):
//     Additional margin around the volume's bounding box to define the search area.
//
// - GDML File Prefix (TString):
//     Prefix of the GDML file to load geometry from.
//
// - Step Fraction (double):
//     Fraction of box dimensions used as step size for scanning the box.
//
// Key Features:
// -------------
// - Box Origin:
//     - Retrieves the origin of the specified volume in local coordinates.
//     - Iterates over the nodes in the top-level geometry to compute the
//       local-to-global transformation matrix for the volume.
//     - Converts the local coordinates to global coordinates in millimeters.
//
// - Box Dimensions:
//     - Defines a 3D box centered on the computed global coordinates.
//     - The box has half-lengths (in mm) derived from the volume's bounding
//       box, extended by the specified margin.
//
// - Volume Search:
//     - Iteratively scans the 3D box, stepping through its dimensions in fine
//       increments starting from the center.
//     - Grid granularity is controlled by the `stepFraction` parameter.
//     - At each grid point, determines the intersecting volume name using the
//       `find_volume_at_point` function.
//
// - Volume Visualization:
//     - Makes intersecting volumes visible in the ROOT OpenGL viewer.
//     - Assigns distinct colors to the nearest volumes for differentiation.
//
// Outputs:
// --------
// - Logs the names of intersecting volumes and their respective intersection points
//   within the box.
// - Launches a ROOT `TBrowser` for interactive exploration of the geometry.
//   (Tip: Use a clipping plane in the viewer to inspect the boundaries.)
//
// How to Run:
// -----------
// Execute the macro in ROOT using the following syntax:
//
// root -l draw_nearby_volumes.C\("[volume_name]", [margin_in_mm], "[gdml_file_prefix]", [step_fraction]\)
//
// Example:
// --------
//
// root -l draw_nearby_volumes.C\("SRetainingWallFoot", 300., "mu2e_40", 0.01\)
//
// Notes:
// ------
// - Adjust the `margin_mm` and `step_fraction` parameters to control
//   the search area and grid resolution.
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


// Utility function to find the volume at a given point
TString find_volume_at_point(double* pos, TGeoManager* geom) {
    // Convert position from mm to cm
    double x = pos[0] / 10.0;
    double y = pos[1] / 10.0;
    double z = pos[2] / 10.0;

    TGeoNode* node = geom->FindNode(x, y, z);
    if (!node) {
        std::cerr << "Position (" << pos[0] << ", " << pos[1] << ", " << pos[2]
                  << ") is outside the world volume." << std::endl;
        return "";
    }

    return clean_volume_names(node->GetName());
}


int draw_nearby_volumes(
    TString targetVolumeName = "SRetainingWallFoot",
    double margin_mm = 300.0,
    TString gdmlFilePrefix = "mu2e_40",
    double stepFraction = 0.05) {

    // *********** IMPORT GEOMETRY *************
    TString gdmlFile = gdmlFilePrefix + ".gdml";
    if (gSystem->AccessPathName(gdmlFile.Data())) {
        std::cerr << "Error: GDML file not found: " << gdmlFile << std::endl;
        return -1;
    }
    TGeoManager* geom = TGeoManager::Import(gdmlFile);

    // ********* FIND THE TARGET VOLUME *********
    TGeoVolume* targetVol = geom->FindVolumeFast(targetVolumeName.Data());
    if (!targetVol) {
        std::cerr << "Error: Volume not found: " << targetVolumeName << std::endl;
        return 1;
    }

    // Hide all volumes initially
    for (const auto obj : *geom->GetListOfVolumes()) {
        TGeoVolume* vol = static_cast<TGeoVolume*>(obj);
        vol->SetVisibility(0);
    }

    // ********* GET BOUNDING BOX PARAMETERS *********
    TGeoBBox* bbox = static_cast<TGeoBBox*>(targetVol->GetShape());
    const double* localBBoxOrigin = bbox->GetOrigin();

    // Convert origin to global coordinates
    double globalBBoxOrigin[3];
    TGeoIterator next(geom->GetTopVolume());
    TGeoNode* current;
    bool originFound = false;
    while ((current = next())) {
        if (current->GetVolume()->GetName() == targetVolumeName) {
            auto globalMatrix = next.GetCurrentMatrix();
            globalMatrix->LocalToMaster(localBBoxOrigin, globalBBoxOrigin);

            // Convert cm to mm
            for (int i = 0; i < 3; ++i) {
                globalBBoxOrigin[i] *= 10.0;
            }
            std::cout << "Volume " << targetVolumeName << " is at global position: ("
                      << globalBBoxOrigin[0] << ", " << globalBBoxOrigin[1] << ", "
                      << globalBBoxOrigin[2] << ") mm" << std::endl;
            originFound = true;
            break;
        }
    }

    if (!originFound) {
        std::cerr << "Error: Unable to find origin of volume " << targetVolumeName
                  << " in global coordinates." << std::endl;
        return -3;
    }

    // Get bounding box dimensions with margin
    double halfLength[3] = {
        bbox->GetDY() * 10.0 + margin_mm, // X in mm
        bbox->GetDZ() * 10.0 + margin_mm, // Y in mm
        bbox->GetDX() * 10.0 + margin_mm  // Z in mm
    };

    // Calculate step size
    int nSteps = static_cast<int>(2.0 / stepFraction);
    double step[3] = {
        2.0 * halfLength[0] / nSteps,
        2.0 * halfLength[1] / nSteps,
        2.0 * halfLength[2] / nSteps
    };

    // ********* FIND NEARBY VOLUMES *********
    std::vector<TString> foundVolumes;
    std::vector<int> colors = {kRed, kOrange, kBlue, kMagenta, kCyan, kSpring};
    int coloredVolumeCount = 0;

    // Loop through grid points
    double gridPoint[3]; // In mm
    for (int ix = 0; ix < nSteps; ++ix) {
        gridPoint[0] = globalBBoxOrigin[0] + ((ix % 2 == 0) ? ix / 2 * step[0] : -(ix + 1) / 2 * step[0]);
        for (int iy = 0; iy < nSteps; ++iy) {
            gridPoint[1] = globalBBoxOrigin[1] + ((iy % 2 == 0) ? iy / 2 * step[1] : -(iy + 1) / 2 * step[1]);
            for (int iz = 0; iz < nSteps; ++iz) {
                gridPoint[2] = globalBBoxOrigin[2] + ((iz % 2 == 0) ? iz / 2 * step[2] : -(iz + 1) / 2 * step[2]);

                TString nearbyVolName = find_volume_at_point(gridPoint, geom);
                if (nearbyVolName != "HallAir" && !nearbyVolName.IsNull()) {
                    if (std::find(foundVolumes.begin(), foundVolumes.end(), nearbyVolName) == foundVolumes.end()) {
                        TGeoVolume* vol = geom->FindVolumeFast(nearbyVolName.Data());
                        if (vol) {
                            foundVolumes.push_back(nearbyVolName);
                            vol->SetVisibility(1);

                            // Assign color to the first few volumes
                            if (coloredVolumeCount < static_cast<int>(colors.size())) {
                                vol->SetLineColor(colors[coloredVolumeCount++]);
                            }
                        }
                        std::cout << "Volume: " << nearbyVolName << " at ("
                                  << gridPoint[0] << ", " << gridPoint[1] << ", "
                                  << gridPoint[2] << ") mm" << std::endl;
                    }
                }
            }
        }
    }

    // Summary of found volumes
    std::cout << foundVolumes.size() << " volumes intersecting the bounding box."
              << std::endl;

    // ********* DRAW VOLUMES *********
    TBrowser* browser = new TBrowser();
    geom->GetVolume("HallAir")->Draw("ogl");

    return 0;
}
