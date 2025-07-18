#ifndef DataProducts_CaloSiPMId_hh
#define DataProducts_CaloSiPMId_hh
//
// Offline identifier of one calorimeter SiPM channel
// Online numbering is in CaloRawSiPMId
//

#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DataProducts/inc/CrystalId.hh"
#include <iosfwd>
#include <bitset>

namespace mu2e {

  class CaloSiPMId{

  public:

    using value_type = CaloConst::CaloSiPMId_type;

    CaloSiPMId():_id(CaloConst::_invalid) {}
    explicit CaloSiPMId(value_type id):_id(id) {}

    value_type id() const { return _id; }
    value_type SiPMLocalId() const { return _id%CaloConst::_nSiPMPerCrystal; }
    CrystalId crystal() const
    { return CrystalId(CrystalId::value_type(_id/CaloConst::_nSiPMPerCrystal)); }
    value_type disk() const;
    bool isValid() const { return _id < CaloConst::_nChannel; }
    bool isCrystal() const { return _id < CaloConst::_nCrystalChannel; }
    bool isPINDiode() const { return _id >= CaloConst::_nCrystalChannel; }
    CaloConst::detType detType() const;

    std::bitset<4> pinDiodeCode() const { return std::bitset<4>(_id-CaloConst::_nCrystalChannel); }
    int pinDiodeDisk() const { return pinDiodeCode()[0]; }
    int pinDiodePhi() const { return pinDiodeCode()[1]; }
    int pinDiodeSphere() const { return pinDiodeCode()[2]; }
    int pinDiodePin() const { return pinDiodeCode()[3]; }

  private:

    value_type _id;

  };

  std::ostream& operator<<(std::ostream& ost,
                           const CaloSiPMId& id );

}
#endif /* DataProducts_CaloSiPMId_hh */
