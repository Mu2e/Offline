#ifndef CaloConditions_CalSimParams_hh
#define CaloConditions_CalSimParams_hh
//
// CalSimParams provides MC calib constants for Readout simulation
//
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DataProducts/inc/CrystalId.hh"
#include "fhiclcpp/ParameterSet.h"
#include <array>
#include <iostream>


namespace mu2e {

  class CalSimParams : virtual public ProditionsEntity {
    public:

      typedef std::shared_ptr<CalSimParams>                        ptr_t;
      typedef std::shared_ptr<const CalSimParams>                  cptr_t;
      typedef std::array<CrystalId,          CaloConst::_nCrystal> IdArray;
      typedef std::array<float,              CaloConst::_nCrystal> floatArray;
      typedef std::array<std::vector<float>, CaloConst::_nCrystal> floatVArray;

      constexpr static const char* cxname = {"CalSimParams"};

      CalSimParams(const IdArray& crystalIds, const floatArray& LRUs,
                   const floatVArray& pePerMeVs,const floatVArray& ADCPerMeVs) :
         ProditionsEntity(cxname),
         _crystalIds(crystalIds),
         _LRUs(LRUs),
         _pePerMeVs(pePerMeVs),
         _ADCPerMeVs(ADCPerMeVs)
      {};

      virtual ~CalSimParams() = default;

      void  print(std::ostream& os) const;
      const std::vector<float> pePerMeVs (const CrystalId& Id) const;
      const std::vector<float> ADCPerMeVs(const CrystalId& Id) const;
      float                    LRU       (const CrystalId& Id) const;

    private:
      IdArray     _crystalIds;
      floatArray  _LRUs;
      floatVArray _pePerMeVs;
      floatVArray _ADCPerMeVs;
  };

}

#endif
