//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// $Id: CosmicTuple.hh,v 1.5 2011/02/03 18:38:27 wasiko Exp $
// $Author: wasiko $
// $Date: 2011/02/03 18:38:27 $
//
// Original author Yury Kolomensky (Rob Kutschke)
//

// Framework includes.
#include <string>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Framework/interface/EDFilter.h"
// Mu2e includes.
#include "ToyDP/inc/StepPointMCCollection.hh"

class TH1F;
class TNtuple;

namespace mu2e {

  class CosmicTuple : public edm::EDFilter {
  public:
    
    explicit CosmicTuple(edm::ParameterSet const& pset);
    virtual ~CosmicTuple() { }

    virtual void beginJob(edm::EventSetup const&);
    virtual bool beginRun(edm::Run &r, edm::EventSetup const& eSetup );
 
    // This is called for each event.
    virtual bool filter(edm::Event& e, edm::EventSetup const&);

  private:

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Cut on the minimum energy.
    double _minimump;
    double _maximump;
    int _minHits;
    int _runNumber;

    // Number of events analyzed.
    int _nAnalyzed;

    // Pointers to histograms & ntuples
    TH1F* _hEventsize;
    TNtuple* _ntupTrk;

  };
  
} // end namespace mu2e
