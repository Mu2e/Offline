//
// Object to perform helix fit to straw hits
//
// $Id: TrkHelixFitIT.hh,v 1.1 2012/05/18 18:14:36 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2012/05/18 18:14:36 $
//
#ifndef TrkHelixFitIT_HH
#define TrkHelixFitIT_HH

#include "BaBar/BaBar.hh"

#include "TrkPatRec/inc/TrkHelixFit.hh"
#include "FastPatternReco/inc/FastPatRecoUtilsAndDataDef.hh"

namespace mu2e 
{
// output struct
  class TrkHelixFitIT:public TrkHelixFit
  {
  public:
// parameter set should be passed in on construction
    explicit TrkHelixFitIT(fhicl::ParameterSet const&);
    virtual ~TrkHelixFitIT();
// main function: given a track definition, find the helix parameters
    bool findHelix(TrkDef const& mytrk,const points3D& pt,TrkHelix& myfit);
  private:
// utlity functions
    void fillXYZP(TrkDef const& mytrk,const points3D& pt, std::vector<XYZP>& xyzp);
    
  };
}
#endif
