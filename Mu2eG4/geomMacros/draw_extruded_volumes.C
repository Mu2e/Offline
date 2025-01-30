//
// Root macro to draw up to 4 extruded volumes with their 3D projections
// and vertices. Rotated extruded solids are also allowed.
//
// Inputs:
// -------
//
//   - The gdml file prefix (omitting .gdml)
//
//   - The number of volumes (from 1 to 4)
//
//   - The volume names (TString)
//     N.B.: The first volume defines the range of 2D projection histograms
//
// Key features:
// -------------
//
//   - Geometry Loading:
//       It imports a GDML geometry file to obtain the required volumes.
//
//   - Volume Verification:
//       It ensures the specified volumes exist and are of type TGeoXtru.
//
//   - Bounding Box Computation:
//       Computes bounding box properties in both local coordinate systems.
//
//   - 2D Projections and Perimeter Visualization:
//        Visualizes the 2D perimeter for each volume in its local coordinates.
//        Provides 2D projections (Top, Side, Front) of the volumes in global
//        coordinates.
//
//   - Interactive 3D Visualization:
//        Opens an interactive browser to explore the specified volumes,
//        with visibility settings customized for clarity.
//
// Output:
// -------
//
//    - Logs volumes vertices.
//
//    - Draw volumes 2D projections
//
//    - Opens a TBrowser for interactive exploration.
//
// How to run:
// -----------
//
//    (from Rob Kutschke's hint (July 2024):
//    To use the OpenGL viewer on Alma9 on mu2e machines you need to:
//       - open a clean shell ( not even mu2einit )
//       - source /cvmfs/scisoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
//       - spack load root/xarsy2v
//
//    root -l draw_extruded_volumes.C\(\"[gdml_file_prefix]\",[number of volumes],\"[1st_volume_name]\",\"[2nd_volume_name]\",\"[3rd_volume_name]\",\"[4th_volume_name]\"\)
//
//    Examples:
//    root -l draw_extruded_volumes.C\(\"mu2e_40\",1,\"SRetainingWallFoot\"\)
//    root -l draw_extruded_volumes.C\(\"mu2e_40\",4,\"SRetainingWallFoot\",\"SRetainingWall\",\"WTSareaCF1notch\",\"WTSareabeam\"\)
//
// Author: S. Di Falco
// Jan 23, 2025
//

#include <iostream>
#include <cmath>
#include "TStyle.h"
#include "TGeoVolume.h"
#include "TGeoXtru.h"
#include "TCanvas.h"
#include "TH2.h"

#define MAXNPOINTS 200
#define MAXNVOLUMES 4
#define nPROJs 3  // number of 2D projections


// Helper function: Load Geometry
TGeoManager* loadGeometry(const TString& gdmlFile) {
  if (gSystem->AccessPathName(gdmlFile.Data())) {
    std::cerr << "Error: GDML file not found: " << gdmlFile << std::endl;
    return nullptr;
  }
  return TGeoManager::Import(gdmlFile);
}

// Helper function: Check if volumes exist
bool validateVolumes(TGeoManager* geom, const TString volNames[], TGeoVolume* vol[], int nVOLUMEs) {
  for (int i = 0; i < nVOLUMEs; i++) {
    vol[i] = geom->FindVolumeFast(volNames[i].Data());
    if (!vol[i]) {
      std::cerr << "Error: Volume not found: " << volNames[i] << std::endl;
      return false;
    }
  }
  return true;
}

// Helper function: Validate if volumes are extruded solids
bool validateExtrudedVolumes(TGeoVolume* vol[], TGeoXtru* xtru[], const TString volNames[], int nVOLUMEs) {
  for (int i = 0; i < nVOLUMEs; i++) {
    xtru[i] = dynamic_cast<TGeoXtru*>(vol[i]->GetShape());
    if (!xtru[i]) {
      std::cerr << "Error: Volume " << volNames[i] << " is not an extruded solid." << std::endl;
      return false;
    }
  }
  return true;
}

// Helper function: Get bounding box information
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
getBoundingBox(TGeoVolume* volume) {
  TGeoBBox* bbox = dynamic_cast<TGeoBBox*>(volume->GetShape());
  if (!bbox) {
    std::cerr << "Error: Bounding box could not be retrieved for volume " << volume->GetName() << std::endl;
    return {{}, {}, {}};
  }

  const double* origin = bbox->GetOrigin();
  std::vector<double> originVec = {origin[0], origin[1], origin[2]};
  std::vector<double> halfLengthVec = {bbox->GetDX(), bbox->GetDY(), bbox->GetDZ()};
  std::vector<double> edgeVec = {
    origin[0] + halfLengthVec[0],
    origin[1] + halfLengthVec[1],
    origin[2] + halfLengthVec[2]};

  return {originVec, halfLengthVec, edgeVec};
}

// Function to project 2D coordinates onto global axes
void projectGlobalCoordinates(const TGeoMatrix* local2global, const TGeoXtru* xtru, std::vector<std::vector<double>>& global_coords, int layer) {
  double local_coo[3];
  double global_coo[3];

  int nVertices=xtru->GetNvert();
  for (int ivert = 0; ivert < nVertices; ivert++) {
    local_coo[0] = xtru->GetX(ivert);
    local_coo[1] = xtru->GetY(ivert);
    local_coo[2] = xtru->GetZ(layer);

    local2global->LocalToMaster(local_coo, global_coo);

    for (int i = 0; i < 3; i++) {
      global_coords[i].push_back(global_coo[i]);
    }
  }
}

// Function to create and draw 2D projections
void draw2DProjections(const std::vector<std::vector<std::vector<double>>>& x_glob, const TString* volNames, vector <int> vol_color, double marker_size, double margin) {

  int nVOLUMEs=x_glob.size();
  TString projName[nPROJs]={"Top","Side","Front"};

  TCanvas* c2D = new TCanvas("c2D", "2D projections", 0, 0, 1800, 600);
  c2D->Divide(nPROJs, 1);

  TString axis_name[3]={"x","y","z"};
  int xaxis[3] = {2, 2, 0}; // Axis mapping
  int yaxis[3] = {0, 1, 1};

  for (int iproj = 0; iproj < nPROJs; iproj++) {
    // Frame initialization
    double xmin = *std::min_element(x_glob[0][xaxis[iproj]].begin(), x_glob[0][xaxis[iproj]].end());
    double xmax = *std::max_element(x_glob[0][xaxis[iproj]].begin(), x_glob[0][xaxis[iproj]].end());
    double ymin = *std::min_element(x_glob[0][yaxis[iproj]].begin(), x_glob[0][yaxis[iproj]].end());
    double ymax = *std::max_element(x_glob[0][yaxis[iproj]].begin(), x_glob[0][yaxis[iproj]].end());

    double maxLength = std::max(xmax-xmin,ymax-ymin);

    double hXmin = xmin - margin;
    double hXmax = xmin + maxLength + margin;
    double hYmin = ymin - margin;
    double hYmax = ymin + maxLength + margin;

    TString hTitle = projName[iproj] + " View; " + axis_name[xaxis[iproj]] + "(cm); " + axis_name[yaxis[iproj]] + "(cm)";
    TH2F* h2proj = new TH2F(Form("h2proj%d", iproj), hTitle.Data(), 5000, hXmin, hXmax, 5000, hYmin, hYmax);
    h2proj->SetStats(0);

    // Graph drawing
    c2D->cd(iproj + 1);
    h2proj->Draw();

    for (int ivol = 0; ivol < nVOLUMEs; ivol++) {
      TGraph* gProj = new TGraph(x_glob[ivol][xaxis[iproj]].size());
      gProj->SetTitle(volNames[ivol]);
      gProj->SetMarkerStyle(20);
      gProj->SetMarkerSize(marker_size);
      gProj->SetLineWidth(2);
      gProj->SetMarkerColor(vol_color[ivol]);
      gProj->SetLineColor(vol_color[ivol]);

      for (size_t i = 0; i < x_glob[ivol][xaxis[iproj]].size(); i++) {
        gProj->SetPoint(i, x_glob[ivol][xaxis[iproj]][i], x_glob[ivol][yaxis[iproj]][i]);
      }

      gProj->Draw("Pl");
    }
  }
}

int draw_extruded_volumes(
                          TString gdmlFilePrefix = "mu2e_40",
                          int nVOLUMEs=4,
                          TString vol1Name="SRetainingWallFoot",
                          TString vol2Name="SRetainingWall",
                          TString vol3Name="WTSareaCF1notch",
                          TString vol4Name="WTSareaCFbeam"
                          ) {

  // Define constants
  gStyle->SetTitleOffset(1.6, "y");
  double margin = 50.; // Margin for the bounding box (in cm)
  std::vector <int> vol_color= {kRed, kBlue, kOrange, kMagenta};
  double marker_size = 1.2;

  // Load geometry
  TString gdmlFile = gdmlFilePrefix + ".gdml";
  TGeoManager* geom = loadGeometry(gdmlFile);
  if (!geom) return -1;

  // Volume names and validation
  TString volNames[MAXNVOLUMES] = {vol1Name, vol2Name, vol3Name, vol4Name};
  TGeoVolume* vol[MAXNVOLUMES];
  if (!validateVolumes(geom, volNames, vol, nVOLUMEs)) return 1;

  // Validate extruded solids
  TGeoXtru* xtru[MAXNVOLUMES];
  if (!validateExtrudedVolumes(vol, xtru, volNames, nVOLUMEs)) return 2;


  // Retrieve bounding box information for each volume
  std::vector<std::vector<double>> local_bbox_origin(nVOLUMEs);
  std::vector<std::vector<double>> local_bbox_halfLength(nVOLUMEs);
  std::vector<std::vector<double>> local_bbox_edge(nVOLUMEs);

  for (int ivol = 0; ivol < nVOLUMEs; ivol++) {
    auto [origin, halfLength, edge] = getBoundingBox(vol[ivol]);

    if (origin.empty() || halfLength.empty() || edge.empty()) {
      std::cerr << "Error: Failed to retrieve bounding box data for volume: " << volNames[ivol] << std::endl;
      return 3;
    }

    local_bbox_origin[ivol] = origin;
    local_bbox_halfLength[ivol] = halfLength;
    local_bbox_edge[ivol] = edge;
  }

  // Draw perimeter and vertices
  TCanvas* cperim[MAXNVOLUMES];
  TGraph* gperim[MAXNVOLUMES];
  TH2F* h2perim[MAXNVOLUMES];
  char cname[100];
  char hname[100];
  char tpoint[100];
  TText *tperim[MAXNVOLUMES][MAXNPOINTS];

  for (int i = 0; i < nVOLUMEs; i++) {
    int nVertices = xtru[i]->GetNvert();
    std::cout << "** " << volNames[i] << " ** has " <<  nVertices << " vertices:" << std::endl;
    std::vector<double> x_loc(nVertices), y_loc(nVertices);

    gperim[i] = new TGraph(nVertices + 1);
    gperim[i]->SetTitle(volNames[i].Data());
    gperim[i]->SetMarkerStyle(20);
    gperim[i]->SetMarkerSize(marker_size);
    gperim[i]->SetMarkerColor(vol_color[i]);
    gperim[i]->SetLineColor(vol_color[i]);
    gperim[i]->SetLineWidth(2);

    for (int j = 0; j < nVertices; j++) {
      x_loc[j] = xtru[i]->GetX(j);
      y_loc[j] = xtru[i]->GetY(j);
      gperim[i]->SetPoint(j, x_loc[j], y_loc[j]);
      std::cout << j << " (" << x_loc[j] << "," << y_loc[j] << ")" << std::endl;
    }
    std::cout << "Bottom z=" << xtru[i]->GetZ(0) << " Top z=" <<  xtru[i]->GetZ(1) << std::endl;
    gperim[i]->SetPoint(nVertices, x_loc[0], y_loc[0]);

    sprintf(cname, "cperim%d", i);
    cperim[i] = new TCanvas(cname, volNames[i].Data(), 600 * i, 0, 500, 500);
    double xmin = *std::min_element(x_loc.begin(), x_loc.end());
    double ymin = *std::min_element(y_loc.begin(), y_loc.end());
    double maxLength = 2. * std::max(local_bbox_halfLength[i][0], local_bbox_halfLength[i][1]);
    sprintf(hname, "h2perim%d", i);
    h2perim[i] = new TH2F(hname, volNames[i] + " perimeter; x(cm); y(cm)",
                          5000, xmin - margin, xmin + maxLength + margin,
                          5000, ymin - margin, ymin + maxLength + margin);
    h2perim[i]->SetStats(0);
    h2perim[i]->Draw();
    gperim[i]->Draw("PL");

    // Add point labels
    for (int ipoint=0;ipoint<nVertices;ipoint++){
      sprintf(tpoint,"%d ",ipoint);
      tperim[i][ipoint]=new TText(x_loc[ipoint],y_loc[ipoint],tpoint);
      tperim[i][ipoint]->SetTextSize(0.04);
      tperim[i][ipoint]->SetTextColor(vol_color[i]);
      tperim[i][ipoint]->Draw();
    }
  }

  // Variables to hold global data
  std::vector<std::vector<std::vector<double>>> x_glob(nVOLUMEs, std::vector<std::vector<double>>(3));

  const TGeoMatrix* local2global;
  TGeoIterator next(geom->GetTopVolume());
  TGeoNode* current;
  int nfound = 0;

  while ((current = next())) {
    TString vol_name = current->GetVolume()->GetName();

    for (int ivol = 0; ivol < nVOLUMEs; ivol++) {
      if (vol_name != volNames[ivol]) continue;

      local2global = next.GetCurrentMatrix();

      // Process bottom layer
      projectGlobalCoordinates(local2global, xtru[ivol], x_glob[ivol], 0);
      // Close perimeter line
      for (int i = 0; i < 3; i++) {
        x_glob[ivol][i].push_back(x_glob[ivol][i][0]);
      }

      // Process top layer
      projectGlobalCoordinates(local2global, xtru[ivol], x_glob[ivol], 1);

      // Get Bounding Box information in Global coordinates
      const double local_origin[3]={local_bbox_origin[ivol][0],
        local_bbox_origin[ivol][1], local_bbox_origin[ivol][2]};
      const double local_edge[3]={local_bbox_edge[ivol][0],
        local_bbox_edge[ivol][1], local_bbox_edge[ivol][2]};

      nfound++;
      break;
    }
    if (nfound == nVOLUMEs) break;
  }

  if (nfound < nVOLUMEs) {
    std::cerr << "Some volumes are missing the coordinate transformation!" << std::endl;
    return 3;
  }

  // Draw projections
  draw2DProjections(x_glob,volNames,vol_color,marker_size,margin);

  // ************* DRAW INTERACTIVE VOLUME ********************

  // Hide all volumes initially
  for (const auto obj : *geom->GetListOfVolumes()) {
    TGeoVolume* vol = (TGeoVolume*)obj;
    vol->SetVisibility(0);
  }

  // Make the wanted volume visible and colored
  for (int ivol=0;ivol<nVOLUMEs;ivol++){
    vol[ivol]->SetVisibility(1);
    vol[ivol]->SetLineColor(vol_color[ivol]);
  }

  TBrowser* b = new TBrowser();
  geom->GetVolume("HallAir")->Draw("ogl");

  return 0;
}
