#ifndef DataProducts_CrystalId_hh
#define DataProducts_CrystalId_hh
//
// Offline identifier of one calorimeter crystal channel
//

#include "Offline/DataProducts/inc/CaloConst.hh"
#include <iosfwd>

namespace mu2e {

  class CrystalId{

  public:

    using value_type=std::uint16_t;

    CrystalId():_id(CaloConst::_invalid) {}
    explicit CrystalId(value_type id):_id(id) {}

    value_type id() const { return _id; }
    value_type disk() const { return _id/CaloConst::_nCrystalPerDisk; };
    bool isValid() const { return _id < CaloConst::_nCrystal; }
    bool isCaphri() const ;
    // this can be used to construct CaloSiPMId
    CaloConst::CaloSiPMId_type SiPMId(CaloConst::SiPMn SiPM01) const
        { return _id*CaloConst::_nSiPMPerCrystal + SiPM01; }

  private:

    value_type _id;

  };

  std::ostream& operator<<(std::ostream& ost,
                           const CrystalId& id );

}
#endif /* DataProducts_CrystalId_hh */
