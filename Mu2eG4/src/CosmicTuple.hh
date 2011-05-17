//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// $Id: CosmicTuple.hh,v 1.6 2011/05/17 15:36:00 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:00 $
//
// Original author Yury Kolomensky (Rob Kutschke)
//

// Framework includes.
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"

#include "art/Framework/Core/EDFilter.h"
// Mu2e includes.
#include "ToyDP/inc/StepPointMCCollection.hh"

class TH1F;
class TNtuple;

namespace mu2e {

  class CosmicTuple : public art::EDFilter {
  public:
    
    explicit CosmicTuple(fhicl::ParameterSet const& pset);
    virtual ~CosmicTuple() { }

    virtual void beginJob(art::EventSetup const&);
    virtual bool beginRun(art::Run &r, art::EventSetup const& eSetup );
 
    // This is called for each event.
    virtual bool filter(art::Event& e, art::EventSetup const&);

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
