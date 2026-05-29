#ifndef CaloConditions_CalSimCrystal_hh
#define CaloConditions_CalSimCrystal_hh
//
// CalSimCrystal provides MC calib constants for Readout simulation
//
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DataProducts/inc/CrystalId.hh"
#include "fhiclcpp/ParameterSet.h"
#include <array>
#include <iostream>


namespace mu2e {

  class CalSimCrystal : virtual public ProditionsEntity {
    public:

      typedef std::shared_ptr<CalSimCrystal>                       ptr_t;
      typedef std::shared_ptr<const CalSimCrystal>                 cptr_t;
      typedef std::array<CrystalId,          CaloConst::_nCrystal> IdArray;
      typedef std::array<float,              CaloConst::_nCrystal> floatArray;
      typedef std::array<std::vector<float>, CaloConst::_nCrystal> floatVArray;

      constexpr static const char* cxname = {"CalSimCrystal"};

      CalSimCrystal(const IdArray& crystalIds, const floatArray& LRUs,
                    const floatVArray& pePerMeVs) :
         ProditionsEntity(cxname),
         _crystalIds(crystalIds),
         _LRUs(LRUs),
         _pePerMeVs(pePerMeVs)
      {};

      virtual ~CalSimCrystal() = default;

      const std::vector<float> pePerMeVs(const CrystalId& Id) const;
      float LRU(const CrystalId& Id) const;
      void  print(std::ostream& os) const;

    private:
      IdArray     _crystalIds;
      floatArray  _LRUs;
      floatVArray _pePerMeVs;
  };

}

#endif
