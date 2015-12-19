//
// base class to resolve hit ambiguities 
//
// $Id: $
// $Author: $
// $Date: $
//
#ifndef WaveformProcess_HH
#define WaveformProcess_HH

#include "RecoDataProducts/inc/RecoCaloDigiCollection.hh"
#include "RecoDataProducts/inc/RecoCaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "MCDataProducts/inc/CaloDigiMC.hh"

#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

namespace mu2e {

  class WaveformProcess {
    public:
// construct from parameter set
#ifndef __GCCXML__
    explicit WaveformProcess(fhicl::ParameterSet const&);
#endif/*__GCCXML__*/
    virtual ~WaveformProcess() = 0;

// resolve a waveform, filling a collecntion of RecoCaloDigi

    virtual void      processWaveform(double ADCToMeV, CaloDigi CaloHit, RecoCaloDigiCollection &RecoCaloHits)  = 0;
  
    virtual void      fillDiag(CaloDigiMC*DigiMC, RecoCaloDigiCollection *RecoCaloHits) = 0;
    
    virtual void      book() = 0;
  };
}

#endif
