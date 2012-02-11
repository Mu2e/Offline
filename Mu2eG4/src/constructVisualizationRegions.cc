// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructVisualizationRegions.hh"

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "G4Color.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Transform3D.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib/exception.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"

#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"

#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"


//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  //================================================================
  void constructVisualizationRegions(const VolumeInfo& worldVolume, const SimpleConfig& config)
  {
    GeomHandle<WorldG4> worldGeom;

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = false; // overlaps are OK
    const bool placePV             = true;

    std::vector<int> visible;
    config.getVectorInt("visregions.boxes.visible", visible, visible);
    if(!visible.empty()) {

      const int nboxes = visible.size();

      std::vector<int> solid; config.getVectorInt("visregions.boxes.solid", solid, nboxes);
      std::vector<std::string> material; config.getVectorString("visregions.boxes.material", material, nboxes);

      std::vector<double> red; config.getVectorDouble("visregions.boxes.color.red", red, nboxes);
      std::vector<double> green; config.getVectorDouble("visregions.boxes.color.green", green, nboxes);
      std::vector<double> blue; config.getVectorDouble("visregions.boxes.color.blue", blue, nboxes);


      std::vector<double> xmin; config.getVectorDouble("visregions.boxes.xmin", xmin, nboxes);
      std::vector<double> ymin; config.getVectorDouble("visregions.boxes.ymin", ymin, nboxes);
      std::vector<double> zmin; config.getVectorDouble("visregions.boxes.zmin", zmin, nboxes);

      std::vector<double> xmax; config.getVectorDouble("visregions.boxes.xmax", xmax, nboxes);
      std::vector<double> ymax; config.getVectorDouble("visregions.boxes.ymax", ymax, nboxes);
      std::vector<double> zmax; config.getVectorDouble("visregions.boxes.zmax", zmax, nboxes);

      for(int ibox = 0; ibox < nboxes; ++ibox) {

        CLHEP::Hep3Vector boxCenterInMu2e((xmax[ibox]+xmin[ibox])/2,
                                          (ymax[ibox]+ymin[ibox])/2,
                                          (zmax[ibox]+zmin[ibox])/2
                                          );

        std::vector<double> boxHalfSize(3);
        boxHalfSize[0] = (xmax[ibox]-xmin[ibox])/2;
        boxHalfSize[1] = (ymax[ibox]-ymin[ibox])/2;
        boxHalfSize[2] = (zmax[ibox]-zmin[ibox])/2;

        std::ostringstream boxname;
        boxname<<"VisualizationBox"<<std::setw(3)<<std::setfill('0')<<ibox;

        nestBox(boxname.str(),
                boxHalfSize,
                findMaterialOrThrow(material[ibox]),
                0,
                boxCenterInMu2e + worldGeom->mu2eOriginInWorld(),
                worldVolume, 0,
                visible[ibox],
                G4Colour(red[ibox], green[ibox], blue[ibox]),
                solid[ibox],
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      }
    }
  }

} // namespace mu2e
