//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// $Id: ReadBack.hh,v 1.1 2010/02/11 15:24:48 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/02/11 15:24:48 $
//
// Original author Rob Kutschke
//

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"

// Mu2e includes.
#include "ToyDP/inc/StepPointMCCollection.hh"

class TH1F;
class TNtuple;

namespace mu2e {

  class Straw;

  class ReadBack : public edm::EDAnalyzer {
  public:
    
    explicit ReadBack(edm::ParameterSet const& pset);
    virtual ~ReadBack() { }

    virtual void beginJob(edm::EventSetup const&);
 
    // This is called for each event.
    void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

    // Cut on the minimum energy.
    double _minimumEnergy;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Number of events analyzed.
    int _nAnalyzed;

    // Pointers to histograms to be filled.
    TH1F* _hRadius;
    TH1F* _hEnergyDep;
    TH1F* _hTime;
    TH1F* _hMultiplicity;
    TH1F* _hDriftDist;
    TH1F* _hxHit;
    TH1F* _hyHit;
    TH1F* _hzHit;
    TH1F* _hHitNeighbours;
    TH1F* _hCheckPointRadius;

    TNtuple* _ntup;

    // Do the work specific to one of the trackers.
    void ReadBack::doLTracker(const edm::Event& event);
    void ReadBack::doITracker(const edm::Event& event);

    // A helper function.
    int countHitNeighbours( Straw const& straw, 
			    edm::Handle<StepPointMCCollection>& hits );

  };
  
} // end namespace mu2e
