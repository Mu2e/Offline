//
// $Id: MCCaloUtilities.hh,v 1.1 2011/01/20 21:28:39 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/01/20 21:28:39 $
//
// Original author Gianni Onorato
//

#ifndef MCCALOUTILITIES_HH
#define MCCALOUTILITIES_HH

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Framework/interface/Event.h"

// Mu2e includes

#include "ToyDP/inc/SimParticleCollection.hh"

namespace mu2e {

  class MCCaloUtilities {

  public:

    MCCaloUtilities();

    ~MCCaloUtilities();

    void setTrackAndRO(const edm::Event & event,
                       SimParticleCollection::key_type track,
                       uint32_t RO);

    void printOutCaloInfo();
    
    bool fromOutside();
    
    bool primary();
    
    bool generated();
    
    int startingVane();

    int getStartingVane(CLHEP::Hep3Vector origin);

    int localVane();

  private:

    uint32_t _localRO;
    uint32_t _localCrystal;
    uint32_t _localVane;
    int _startingVane;
    bool _fromOutside, _primary, _generated;

  };

}

#endif
