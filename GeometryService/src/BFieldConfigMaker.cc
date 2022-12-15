#include "Offline/GeometryService/inc/BFieldConfigMaker.hh"

#include <iostream>
#include <set>
#include <string>

#include "cetlib_except/exception.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

#include "Offline/BFieldGeom/inc/BFieldConfig.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

    namespace {
        void noLegacyFormatMapSpecs(const SimpleConfig& config) {
            vector<string> keysToCheck;
            keysToCheck.push_back("bfield.files");
            keysToCheck.push_back("bfield.dsFile");
            keysToCheck.push_back("bfield.psFile");
            keysToCheck.push_back("bfield.tsuFile");
            keysToCheck.push_back("bfield.tsdFile");

            for (vector<string>::const_iterator i = keysToCheck.begin(); i != keysToCheck.end();
                 ++i) {
                if (config.hasName(*i)) {
                    throw cet::exception("GEOM") << "Found obsolete config file parameter " << *i
                                                 << ".  Please use bfield.innerMaps and "
                                                    "bfield.outerMaps to specify G4BL maps."
                                                 << "Maps are not loaded.\n";
                }
            }
        }
    }  // namespace

    BFieldConfigMaker::BFieldConfigMaker(const SimpleConfig& config, const Beamline& beamg)
        : bfconf_(new BFieldConfig()) {
        bfconf_->writeBinaries_ = config.getBool("bfield.writeG4BLBinaries", false);
        bfconf_->verbosityLevel_ = config.getInt("bfield.verbosityLevel");
        bfconf_->flipBFieldMaps_ = config.getBool("bfield.flipMaps", false);

        bfconf_->scaleFactor_ = config.getDouble("bfield.scaleFactor", 1.0);

        bfconf_->dsFieldForm_ = BFieldConfig::DSFieldModel(
            config.getInt("detSolFieldForm", BFieldConfig::dsModelUniform));

        const double bz = config.getDouble("toyDS.bz", 0.);
        bfconf_->dsUniformValue_ = CLHEP::Hep3Vector(0., 0., bz * bfconf_->scaleFactor_);

        const double grad = config.getDouble("toyDS.gradient", 0.);
        bfconf_->dsGradientValue_ = CLHEP::Hep3Vector(0., 0., grad * bfconf_->scaleFactor_);

        BFInterpolationStyle style(config.getString("bfield.interpolationStyle", "trilinear"));
        bfconf_->interpStyle_ = style;

        const string format = config.getString("bfield.format","G4BL");

        // Load in the optional formatList, if you're mixing parametric and G4BL
        vector<string> mapTypeList;
        config.getVectorString("bfield.formatList", mapTypeList, vector<string>());
        if (!mapTypeList.empty()) {
            for (const string& mt : mapTypeList) {
                BFMapType bfmt;
                if (mt == "G4BL") {
                    bfmt = BFMapType::G4BL;
                    bfconf_->mapTypeList_.push_back(bfmt);
                } else if (mt == "PARAM") {
                    bfmt = BFMapType::PARAM;
                    bfconf_->mapTypeList_.push_back(bfmt);
                } else {
                    throw cet::exception("GEOM")
                        << "The bfield.formatList must specify a vector of either `G4BL` and/or "
                           "`PARAM`. "
                        << "The specified format, `" << mt << "`, is not valid."
                        << "\n";
                }
            }
        }

        if (format == "G4BL") {
            // These maps require torus radius of 2929 mm
            const double requiredTorusRadius = 2929.0 * CLHEP::mm;
            if (fabs(beamg.getTS().torusRadius() - requiredTorusRadius) > 0.1) {
                throw cet::exception("GEOM")
                    << "The G4BL magnetic field files require torus radius of 2929 mm."
                    << " The value used by geometry is " << beamg.getTS().torusRadius()
                    << " Maps are not loaded.\n";
            }
            bfconf_->mapType_ = BFMapType::G4BL;

            // Make sure we don't have any leftover bits from old geometry format,
            // to avoid user confusion.
            noLegacyFormatMapSpecs(config);

            config.getVectorString("bfield.innerMaps", bfconf_->innerMapFiles_);
            config.getVectorString("bfield.outerMaps", bfconf_->outerMapFiles_);

        } else if (format == "PARAM") {
            // These maps require torus radius of 2929 mm
            const double requiredTorusRadius = 2929.0 * CLHEP::mm;
            if (fabs(beamg.getTS().torusRadius() - requiredTorusRadius) > 0.1) {
                throw cet::exception("GEOM")
                    << "The G4BL magnetic field files require torus radius of 2929 mm."
                    << " The value used by geometry is " << beamg.getTS().torusRadius()
                    << " Maps are not loaded.\n";
            }
            bfconf_->mapType_ = BFMapType::PARAM;

            // Make sure we don't have any leftover bits from old geometry format,
            // to avoid user confusion.
            noLegacyFormatMapSpecs(config);

            config.getVectorString("bfield.innerMaps", bfconf_->innerMapFiles_);
            config.getVectorString("bfield.outerMaps", bfconf_->outerMapFiles_);

        } else {
            throw cet::exception("GEOM")
                << "Unknown format of file with magnetic field maps: " << format << "\n";
        }

        if (!bfconf_->mapTypeList_.empty() &&
            bfconf_->mapTypeList_.size() !=
                bfconf_->innerMapFiles_.size() + bfconf_->outerMapFiles_.size()) {
            throw cet::exception("GEOM") << "If formatList is specified, it must have number of "
                                            "entries equal to number of maps. "
                                         << "\n";
        }
    }


}  // namespace mu2e
