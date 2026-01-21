#ifndef RecoDataProducts_CaloTrigSeed_hh
#define RecoDataProducts_CaloTrigSeed_hh

#include <iostream>
#include <ostream>
namespace mu2e {

  class CaloTrigSeed
  {
      public:

        CaloTrigSeed(): _crystalId(0),_epeak(0.),_tpeak(0.),_rpeak(0.),_ring1max(0.),_ring1max2(0.),_ring2max(0.),_cluenergy(0.),_clutime(0.),_clucogx(0.),_clucogy(0.)
        {}

        CaloTrigSeed(unsigned int crystalId, float epeak, float tpeak, float rpeak, float ring1max, float ring1max2, float ring2max,
                     float cluenergy, float clutime, float clucogx, float clucogy)  :
           _crystalId(crystalId),_epeak(epeak),_tpeak(tpeak),_rpeak(rpeak),_ring1max(ring1max),_ring1max2(ring1max2),_ring2max(ring2max),
           _cluenergy(cluenergy),_clutime(clutime),_clucogx(clucogx),_clucogy(clucogy)
        {         }

        void print(std::ostream& ost = std::cout) const;

        unsigned  crystalid()  const { return _crystalId; }
        float     epeak()      const { return _epeak;}
        float     tpeak()      const { return _tpeak;}
        float     rpeak()      const { return _rpeak;}
        float     ring1max()   const { return _ring1max;}
        float     ring1max2()  const { return _ring1max2;}
        float     ring2max()   const { return _ring2max;}
        float     cluenergy()  const { return _cluenergy;}
        float     clutime()    const { return _clutime;}
        float     clucogx()    const { return _clucogx;}
        float     clucogy()    const { return _clucogy;}


     private:
        unsigned   _crystalId;
        float      _epeak;
        float      _tpeak;
        float      _rpeak;
        float      _ring1max;
        float      _ring1max2;
        float      _ring2max;
        float      _cluenergy;
        float      _clutime;
        float      _clucogx;
        float      _clucogy;
  };


  using CaloTrigSeedCollection = std::vector<mu2e::CaloTrigSeed>;
  using CaloTrigSeedPtrCollection = std::vector<art::Ptr<mu2e::CaloTrigSeed>>;
}

#endif
