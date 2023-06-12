#ifndef CaloConditions_CaloDAQMap_hh
#define CaloConditions_CaloDAQMap_hh

//
// CaloDAQMap converts between offline and readout SiPM channel numbers
//

#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/DataProducts/inc/CaloRawSiPMId.hh"
#include "fhiclcpp/ParameterSet.h"
#include <array>
#include <iostream>

namespace mu2e {

  class CaloDAQMap : virtual public ProditionsEntity {
  public:

    typedef std::shared_ptr<CaloDAQMap> ptr_t;
    typedef std::shared_ptr<const CaloDAQMap> cptr_t;

    typedef std::array<CaloSiPMId,CaloConst::_nRawChannel> RawArray;
    typedef std::array<CaloRawSiPMId,CaloConst::_nChannel> OfflineArray;

    //CaloDAQMap():_name("CaloDAQMap") {}
    constexpr static const char* cxname = {"CaloDAQMap"};

    // construct with constants, then some values are computed and filled below
    CaloDAQMap(const RawArray& raw2Offline, const OfflineArray& offline2Raw) :
       ProditionsEntity(cxname),
       _raw2Offline(raw2Offline),_offline2Raw(offline2Raw) {};

    virtual ~CaloDAQMap() {}

    CaloSiPMId offlineId(CaloRawSiPMId rawId) const;
    CaloRawSiPMId rawId(CaloSiPMId offId) const;

    void print(std::ostream& os) const;

  private:

    RawArray _raw2Offline;
    OfflineArray _offline2Raw;

  };

}

#endif

