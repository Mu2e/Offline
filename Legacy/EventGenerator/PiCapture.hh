#ifndef EventGenerator_PiCapture_hh
#define EventGenerator_PiCapture_hh
//
//
// Generate photons from pi- capture on Al nuclei.
// Based on Ivano Sarra's work described in Mu2e doc 665-v2
//
// $Id: PiCapture.hh,v 1.27 2013/05/31 18:07:29 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 18:07:29 $
//
// Original author Rob Kutschke, P. Shanahan
//

//C++ includes
#include<memory>

// Mu2e includes
#include "EventGenerator/inc/FoilParticleGenerator.hh"
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// CLHEP includes
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "EventGenerator/inc/PiCaptureEffects.hh"

// Forward declarations outside of namespace mu2e.
class TH1D;
namespace art{
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;
  class DetectorSystem;

  class PiCapture: public GeneratorBase{

  public:
    PiCapture( art::Run& run, const SimpleConfig& config );
    virtual ~PiCapture();

    virtual void generate( GenParticleCollection&  );

  private:

    // Start: parameters from run time configuration/

    double _mean;         //< mean per event
    double _elow;         //< lower photon energy
    double _ehi;          //< upper photon energy
    double _probInternalConversion;
    bool _PStoDSDelay;
    bool _pPulseDelay;
    double _pPulseShift;

    // Activate the folding procedure on generation time. Default is on
    bool _timeFolding;

    // Select the position, type and time type for the generation
    std::string _foilGen;
    std::string _posGen;
    std::string _timeGen;

    int    _nbins;        //< number of bins in photon energy pdf
    bool   _doHistograms; // Enable/disable histograms.
    // End: parameters from run time configuration/



    //Class object to pick up random position and time
    std::unique_ptr<FoilParticleGenerator> _fGenerator;


    // Random number distributions
    CLHEP::RandPoissonQ _randPoissonQ;

    std::string _STfname;

    const DetectorSystem *_detSys;

    // Histograms.
    TH1D* _hMultiplicity;
    TH1D* _hEPhot;
    TH1D* _hEPhotZ;
    TH1D* _hEElect;
    TH1D* _hEElectZ;
    TH1D* _hInternalFraction;
    TH1D* _hzPos;
    TH1D* _hcz;
    TH1D* _hphi;
    TH1D* _ht;
    TH1D* _hmudelay;
    TH1D* _hpulsedelay;
    TH1D* _hFoilNumber;

    PiCaptureEffects* _piCaptureInfo;

  };

} // end namespace mu2e,

#endif /* EventGenerator_PiCapture_hh */
