// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructVisualizationRegions.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Geant4/G4Box.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Transform3D.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "cetlib_except/exception.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"

#include "Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"

#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"

#define AGDEBUG(stuff) \
    std::cerr << "AG: " << __FILE__ << ", line " << __LINE__ << ": " << stuff << std::endl;
//#define AGDEBUG(stuff)

namespace mu2e {

    //================================================================
    // reuse same code for both inner and outer maps
    void processFieldMaps(const VolumeInfo& worldVolume,
                          const SimpleConfig& config,
                          const std::string& configPrefix,
                          const BFieldManager::MapContainerType& maps) {

        const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
        geomOptions->loadEntry( config, configPrefix + "visible", configPrefix);

        const bool isVisible            = geomOptions->isVisible(configPrefix + "visible"); 
        const bool isSolid              = geomOptions->isSolid(configPrefix + "visible"); 
        const bool forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible(configPrefix + "visible"); 
        const bool doSurfaceCheck       = geomOptions->doSurfaceCheck(configPrefix + "visible");
        const bool placePV              = geomOptions->placePV(configPrefix + "visible");
        const int  verbosityLevel       = config.getInt("visregions.verbosityLevel", 0);

        if (isVisible) {
            GeomHandle<WorldG4> worldGeom;

            std::vector<double> red;
            config.getVectorDouble(configPrefix + ".color.red", red);
            if (red.empty()) {
                throw cet::exception("GEOM")
                    << "constructVisualizationRegions(): " << configPrefix << "color.red is empty";
            }
            std::vector<double> green;
            config.getVectorDouble(configPrefix + ".color.green", green, red.size());
            std::vector<double> blue;
            config.getVectorDouble(configPrefix + ".color.blue", blue, red.size());

            G4Material* material = findMaterialOrThrow(config.getString(configPrefix + "material"));

            std::vector<double> squashY;
            config.getVectorDouble(configPrefix + ".squashY", squashY, squashY);

            int mapNumber(0);
            // This is a temp hack and will only work for grid-like maps
            for (BFieldManager::MapContainerType::const_iterator i = maps.begin(); i != maps.end();
                 ++i, ++mapNumber) {
                const BFGridMap& field = dynamic_cast<const BFGridMap&>(**i);

                CLHEP::Hep3Vector mapCenterInMu2e(
                    field.xmin() + field.dx() * (field.nx() - 1) / 2,

                    squashY.empty() ? field.ymin() + field.dy() * (field.ny() - 1) / 2
                                    : (squashY[1] + squashY[0]) / 2.,

                    field.zmin() + field.dz() * (field.nz() - 1) / 2);

                std::vector<double> mapHalfSize(3);
                mapHalfSize[0] = field.dx() * (field.nx() - 1) / 2.;
                mapHalfSize[1] = squashY.empty() ? field.dy() * (field.ny() - 1) / 2.
                                                 : (squashY[1] - squashY[0]) / 2.;
                mapHalfSize[2] = field.dz() * (field.nz() - 1) / 2.;

                int icolor = mapNumber % red.size();
                G4Colour color(red[icolor], green[icolor], blue[icolor]);

                std::ostringstream boxname;
                boxname << "VisualizationMapBox_" << field.getKey();

                if (verbosityLevel) {
                    std::cout << __func__ << ": adding " << boxname.str() << " = {";
                    std::copy(mapHalfSize.begin(), mapHalfSize.end(),
                              std::ostream_iterator<double>(std::cout, ", "));
                    std::cout << "} at " << mapCenterInMu2e << std::endl;
                }

                nestBox(boxname.str(), mapHalfSize, material, 0,
                        mapCenterInMu2e + worldGeom->mu2eOriginInWorld(), worldVolume, 0,
                        isVisible, color, isSolid, forceAuxEdgeVisible, placePV,
                        doSurfaceCheck);
            }
        }
    }

    //================================================================
    void constructVisualizationRegions(const VolumeInfo& worldVolume, const SimpleConfig& config) {
        //----------------------------------------------------------------
        GeomHandle<BFieldManager> fieldMgr;
        processFieldMaps(worldVolume, config, "visregions.innerFieldMaps",
                         fieldMgr->getInnerMaps());
        processFieldMaps(worldVolume, config, "visregions.outerFieldMaps",
                         fieldMgr->getOuterMaps());

        //----------------------------------------------------------------
        // Add the arbitrary boxes
        const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
        geomOptions->loadEntry( config, "visregionsBoxes", "visregions.boxes");

        const bool isVisible            = geomOptions->isVisible("visregionsBoxes"); 
        const bool isSolid              = geomOptions->isSolid("visregionsBoxes"); 
        const bool forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("visregionsBoxes"); 
        const bool doSurfaceCheck       = geomOptions->doSurfaceCheck("visregionsBoxes");
        const bool placePV              = geomOptions->placePV("visregionsBoxes");
        const int  verbosityLevel       = config.getInt("visregions.verbosityLevel", 0);

        if (isVisible) {
            GeomHandle<WorldG4> worldGeom;

            G4Material* material =
                findMaterialOrThrow(config.getString("visregions.boxes.material"));

            std::vector<double> red;
            config.getVectorDouble("visregions.boxes.color.red", red);
            if (red.empty()) {
                throw cet::exception("GEOM")
                    << "constructVisualizationRegions(): visregions.boxes.color.red is empty";
            }
            std::vector<double> green;
            config.getVectorDouble("visregions.boxes.color.green", green, red.size());
            std::vector<double> blue;
            config.getVectorDouble("visregions.boxes.color.blue", blue, red.size());

            std::vector<double> xmin;
            config.getVectorDouble("visregions.boxes.xmin", xmin);
            const int nboxes = xmin.size();
            std::vector<double> ymin;
            config.getVectorDouble("visregions.boxes.ymin", ymin, nboxes);
            std::vector<double> zmin;
            config.getVectorDouble("visregions.boxes.zmin", zmin, nboxes);

            std::vector<double> xmax;
            config.getVectorDouble("visregions.boxes.xmax", xmax, nboxes);
            std::vector<double> ymax;
            config.getVectorDouble("visregions.boxes.ymax", ymax, nboxes);
            std::vector<double> zmax;
            config.getVectorDouble("visregions.boxes.zmax", zmax, nboxes);

            for (int ibox = 0; ibox < nboxes; ++ibox) {
                CLHEP::Hep3Vector boxCenterInMu2e((xmax[ibox] + xmin[ibox]) / 2,
                                                  (ymax[ibox] + ymin[ibox]) / 2,
                                                  (zmax[ibox] + zmin[ibox]) / 2);

                std::vector<double> boxHalfSize(3);
                boxHalfSize[0] = (xmax[ibox] - xmin[ibox]) / 2;
                boxHalfSize[1] = (ymax[ibox] - ymin[ibox]) / 2;
                boxHalfSize[2] = (zmax[ibox] - zmin[ibox]) / 2;

                int icolor = ibox % red.size();
                G4Colour color(red[icolor], green[icolor], blue[icolor]);

                std::ostringstream boxname;
                boxname << "VisualizationBox" << std::setw(3) << std::setfill('0') << ibox;

                if (verbosityLevel) {
                    std::cout << __func__ << ": adding " << boxname.str() << " = {";
                    std::copy(boxHalfSize.begin(), boxHalfSize.end(),
                              std::ostream_iterator<double>(std::cout, ", "));
                    std::cout << "} at " << boxCenterInMu2e << std::endl;
                }

                nestBox(boxname.str(), boxHalfSize, material, 0,
                        boxCenterInMu2e + worldGeom->mu2eOriginInWorld(), worldVolume, 0,
                        isVisible, color, isSolid, forceAuxEdgeVisible, placePV,
                        doSurfaceCheck);
            }
        }
    }

}  // namespace mu2e
