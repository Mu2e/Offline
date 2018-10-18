//
// Class to interface with track quality estimator TrkQual
// Original author: Dave Brown LBNL Feb 2017
//
#include "TrkAdapter/inc/TrkQualInter.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include <string>
using std::string;
namespace mu2e {
  TrkQualInter::TrkQualInter(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0))
  {
// initialize TrkQual MVA
    fhicl::ParameterSet mvapset = pset.get<fhicl::ParameterSet>("TrkQualMVA",fhicl::ParameterSet());
    mvapset.put<string>("MVAWeights",pset.get<string>("TrkQualWeights","TrkDiag/test/TrkQual.weights.xml"));
    _trkqualmva.reset(new MVATools(mvapset));
    _trkqualmva->initMVA();
    if(_debug>0)_trkqualmva->showMVA();
  }


  void TrkQualInter::computeMVA(TrkQual& trkqual) const {
    trkqual.setMVAValue(_trkqualmva->evalMVA(trkqual.values()));
  }

  void TrkQualInter::fillInput(KalSeed const& kseed, TrkQual& trkqual) const {
    
  }

}

