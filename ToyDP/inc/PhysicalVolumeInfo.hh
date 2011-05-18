#ifndef ToyDP_PhysicalVolumeInfo_hh
#define ToyDP_PhysicalVolumeInfo_hh

//
// Persistable information about a G4 Physical Volume.
//
// $Id: PhysicalVolumeInfo.hh,v 1.4 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke
//
//

#include <string>

namespace mu2e {

  struct PhysicalVolumeInfo {

    // This c'tor is required for ROOT.
    PhysicalVolumeInfo(){};

    PhysicalVolumeInfo( const std::string& pname,
                        uint32_t           pcopyNo ):
      _name(pname),
      _copyNo(pcopyNo){}

    // Accept compiler generated versions of the
    // destructor, copy constructor and the assignment
    // operator.

    // Accessors

    std::string const& name()   const { return _name;  }
    uint32_t           copyNo() const { return _copyNo;}

  private:
    std::string _name;
    uint32_t _copyNo;

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const PhysicalVolumeInfo& vol ){
    ost << "( "
        << vol.name() << ", "
        << vol.copyNo()
        << " )";
    return ost;
  }


}

#endif /* ToyDP_PhysicalVolumeInfo_hh */
