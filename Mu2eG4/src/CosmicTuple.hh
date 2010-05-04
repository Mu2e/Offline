//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// $Id: CosmicTuple.hh,v 1.1 2010/05/04 01:33:46 yury Exp $
// $Author: yury $
// $Date: 2010/05/04 01:33:46 $
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

// Mu2e includes.
#include "ToyDP/inc/StepPointMCCollection.hh"

using namespace std;

class TH1F;
class TNtuple;

namespace mu2e {

  class CosmicTuple : public edm::EDAnalyzer {
  public:
    
    explicit CosmicTuple(edm::ParameterSet const& pset);
    virtual ~CosmicTuple() { }

    virtual void beginJob(edm::EventSetup const&);
 
    // This is called for each event.
    void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Cut on the minimum energy.
    double _minimumEnergy;

    // Number of events analyzed.
    int _nAnalyzed;

    // Pointers to histograms to be filled.
    TNtuple* _ntupTrk;

  };
  
} // end namespace mu2e
