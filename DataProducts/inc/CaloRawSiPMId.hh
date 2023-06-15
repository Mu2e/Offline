#ifndef DataProducts_CaloRawSiPMId_hh
#define DataProducts_CaloRawSiPMId_hh
//
//
// Online or raw data identifier of one calorimeter SiPM channel
// Offline number identifier is in CaloSiPMId
//

#include "Offline/DataProducts/inc/CaloConst.hh"
#include <iosfwd>

namespace mu2e {

  class CaloRawSiPMId{

  public:

    using value_type=std::uint16_t;

    CaloRawSiPMId():_id(CaloConst::_invalid) {}
    explicit CaloRawSiPMId(value_type id):_id(id) {}

    value_type id() const { return _id; }
    value_type dirac() const { return _id/CaloConst::_nChPerDIRAC; }
    value_type ROCchannel() const { return _id%CaloConst::_nChPerDIRAC; }

    bool isValid() const { return _id < CaloConst::_nRawChannel; }

  private:

    value_type _id;

  };

  std::ostream& operator<<(std::ostream& ost,
                           const CaloRawSiPMId& id );

}
#endif /* DataProducts_CaloRawSiPMId_hh */
