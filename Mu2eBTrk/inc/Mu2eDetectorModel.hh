//--------------------------------------------------------------------------
// Name:
//   Mu2eDetectorModel: top-level detector model for Mu2e 
//
//        Copyright (C) 2015    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 23 Dec 2015
//------------------------------------------------------------------------
#ifndef Mu2eDetectorModel_hh
#define Mu2eDetectorModel_hh
// mu2e includes
#include "fhiclcpp/ParameterSet.h"
// BTrk includes
// Mu2eBTrk includes
#include "DataProducts/inc/StrawId.hh"
#include "Mu2eBTrk/inc/DetStrawElem.hh"
#include "Mu2eBTrk/inc/DetStrawType.hh"
#include "Mu2eInterfaces/inc/Detector.hh"
// c++ includes
#include <map>
#include <string>

namespace mu2e {
  class TTracker;
  
  class Mu2eDetectorModel : public Detector  {
    public:
// construct with parameter set.  Ignored for now
      Mu2eDetectorModel(fhicl::ParameterSet const& pset, TTracker const& tt);
      virtual ~Mu2eDetectorModel();
// given a straw, or index, find the associated element.  This will throw
// if no match is found
    const DetStrawElem* strawElem(Straw const& straw) const;
    const DetStrawElem* strawElem(StrawId const& strawid) const;
    private:
      // map between straw index and detector elements 
      std::map<StrawId,DetStrawElem*> _strawmap;
      // types for straw elements
      DetStrawType* _strawtype; // straw materials description
      // materials
      std::string _gasmatname, _wallmatname, _wirematname;
  };
}

#endif
