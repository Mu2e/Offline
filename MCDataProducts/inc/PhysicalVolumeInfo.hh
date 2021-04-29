#ifndef MCDataProducts_PhysicalVolumeInfo_hh
#define MCDataProducts_PhysicalVolumeInfo_hh

//
// Persistable information about a G4 Physical Volume.
//
//
// Original author Rob Kutschke
//
//

#include <iostream>
#include <string>

namespace mu2e {

  struct PhysicalVolumeInfo {

    // This c'tor is required for ROOT.
    PhysicalVolumeInfo(){}

    PhysicalVolumeInfo( const std::string& pname,
                        unsigned           pcopyNo,
                        const std::string& materialName)
      : _name(pname)
      , _copyNo(pcopyNo)
      , _materialName(materialName)
    {}

    // Accept compiler generated versions of the
    // destructor, copy constructor and the assignment
    // operator.

    // Accessors

    std::string const& name()   const { return _name;  }
    unsigned           copyNo() const { return _copyNo;}
    std::string const& materialName() const { return _materialName;  }

  private:
    std::string _name;
    unsigned _copyNo=0;
    std::string _materialName;
  };

  // Shift left (printing) operator.
  std::ostream& operator<<(std::ostream& ost, const PhysicalVolumeInfo& vol);

  bool operator==(const PhysicalVolumeInfo& a, const PhysicalVolumeInfo& b);
  inline bool operator!=(const PhysicalVolumeInfo& a, const PhysicalVolumeInfo& b) { return !(a==b); }
}

#endif /* MCDataProducts_PhysicalVolumeInfo_hh */
