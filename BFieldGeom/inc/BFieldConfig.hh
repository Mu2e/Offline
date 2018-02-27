//
// Complete configuration of Mu2e magnetic fields as specified by user.
//
// Andrei Gaponenko, 2012, cut-and-pasting code from BFieldManager.

#ifndef BFieldGeom_BFieldConfig_hh
#define BFieldGeom_BFieldConfig_hh

#include <string>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

#include "BFieldGeom/inc/BFInterpolationStyle.hh"
#include "BFieldGeom/inc/BFMapType.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

    class BFieldConfig : virtual public Detector {
       public:
        // The constructor is private, an instance should
        // be obtained from the maker
        friend class BFieldConfigMaker;

        // G4BL vs GMC vs PARAM maps
        BFMapType mapType() const { return mapType_; }
        std::vector<BFMapType> mapTypeList() const { return mapTypeList_; }

        // Linear vs quadratic interpolation.
        BFInterpolationStyle interpolationStyle() const { return interpStyle_; }

        // The order of entries is significant for outer maps, but
        // irrelevant for inner ones.  For simplicity we'll use the same
        // container type to represent both cases.
        // Entries are unique over the union of inner and outer maps.
        typedef std::vector<std::string> FileSequenceType;
        const FileSequenceType& innerMapFiles() const { return innerMapFiles_; }
        const FileSequenceType& outerMapFiles() const { return outerMapFiles_; }

        typedef std::vector<std::vector<int> > GMCDimSequence;
        const GMCDimSequence& gmcDimensions() const { return gmcDimensions_; }

        double scaleFactor() const { return scaleFactor_; }

        // Models of the DS magnetic field:
        // 0 - whole DS uses the field map.
        // 1 - upstream uses the full field map; downstream uses a uniform field.
        // 2 - whole DS uses a uniform field.
        enum DSFieldModel { dsModelFull, dsModelSplit, dsModelUniform };

        // Special field for the detector solenoid to be used instead of the maps
        DSFieldModel dsFieldForm() const { return dsFieldForm_; }

        // The uniform field component in the DS for the special case
        const CLHEP::Hep3Vector& getDSUniformValue() const { return dsUniformValue_; }

        // The gradient field in the DS for the special case
        bool useDSGradientValue() const;

        const CLHEP::Hep3Vector& getDSGradientValue() const { return dsGradientValue_; }

        // to trigger the map-writing hack inside the BFieldManagerMaker code.
        bool writeBinaries() const { return writeBinaries_; }

        int verbosityLevel() const { return verbosityLevel_; }

        bool flipBFieldMaps() const { return flipBFieldMaps_; }

       private:
        BFieldConfig()
            : scaleFactor_(1.), writeBinaries_(false), verbosityLevel_(1), flipBFieldMaps_(false) {}

        // GMC, G4BL or possible future types.
        BFMapType mapType_;
        std::vector<BFMapType> mapTypeList_;

        // Linear, MECO style quadratic or possible future types.
        BFInterpolationStyle interpStyle_;

        FileSequenceType innerMapFiles_;
        FileSequenceType outerMapFiles_;

        // Extra info is needed to interpred GMC maps
        GMCDimSequence gmcDimensions_;

        double scaleFactor_;

        DSFieldModel dsFieldForm_;

        // Special case: uniform field in the DS.
        CLHEP::Hep3Vector dsUniformValue_;

        // Special case: gradient field in the DS
        CLHEP::Hep3Vector dsGradientValue_;

        bool writeBinaries_;
        int verbosityLevel_;
        bool flipBFieldMaps_;
    };

}  // namespace mu2e

#endif /* BFieldGeom_BFieldConfig_hh */
