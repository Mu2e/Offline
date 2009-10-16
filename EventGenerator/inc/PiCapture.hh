#ifndef PiCapture_HH
#define PiCapture_HH

//
//
// Based on Ivano Sarra's work described in Mu2e doc 665-v2
// 
// $Id: PiCapture.hh,v 1.2 2009/10/16 04:20:52 shanahan Exp $
// $Author: shanahan $ 
// $Date: 2009/10/16 04:20:52 $
//
// Original author Rob Kutschke, P. Shanahan
// 


#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "CLHEP/Random/RandGeneral.h"

using CLHEP::RandGeneral;

class TH1D;

namespace edm{
  class Run;
}


namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class PiCapture: public GeneratorBase{

  public:
    PiCapture( edm::Run& run, const SimpleConfig& config );
    virtual ~PiCapture();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    const double EPhotFunc(const double x);

    TH1D* _piCaptureMultiplicity;
    TH1D* _piCaptureEPhot;
    TH1D* _piCaptureEPhotZ;

    RandomUnitSphere _randomUnitSphere;
    RandGeneral* _funcGen;

    double _mean; //< mean per event
    double _elow; //< lower photon energy 
    double _ehi; //< upper photon energy 
    double _bindE; //< energy bin width for generator
    int _nbins; //< number of bins in photon energy pdf

  };

} // end namespace mu2e,

#endif


