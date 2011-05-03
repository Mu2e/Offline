//
// An EDAnalyzer module that reads back the hits created by G4 and makes 
// histograms, ntuples and TGraphs.
//
// $Id: ReadBack.hh,v 1.19 2011/05/03 04:28:25 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/03 04:28:25 $
//
// Original author Rob Kutschke
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

class TH1F;
class TH2F;
class TGraph;
class TNtuple;

namespace mu2e {

  class Straw;

  class ReadBack : public edm::EDAnalyzer {
  public:
    
    explicit ReadBack(edm::ParameterSet const& pset);
    virtual ~ReadBack() { }

    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();
 
    // This is called for each event.
    virtual void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

    // Start: run time parameters

    // Diagnostics printout level
    int _diagLevel;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Module label of the generator module that was passed as input to G4.
    std::string _generatorModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Module which made the CaloCrystalHits
    std::string _caloCrystalHitsMaker;
 
    // Name of the stopping target StepPoint collection
    std::string _targetStepPoints;

    // Name of the CRSScintillatorBar(CRV) StepPoint collection
    std::string _crvStepPoints;

    // Cut on the minimum energy.
    double _minimumEnergy;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Limit the size of the TGraph.
    int _xyHitsMax;

    // End: run time parameters

    // Number of events analyzed.
    int _nAnalyzed;

    // Pointers to histograms, ntuples, TGraphs.
    TH1F* _hRadius;
    TH1F* _hEnergyDep;
    TH1F* _hTime;
    TH1F* _hMultiplicity;
    TH1F* _hDriftDist;
    TH1F* _hDriftDistW;
    TH1F* _hxHit;
    TH1F* _hyHit;
    TH1F* _hzHit;
    TH1F* _hHitNeighbours;
    TH1F* _hCheckPointRadius;
    TH1F* _hCheckPointRadiusW;
    TH1F* _hCheckPointWireZ;
    TH1F* _hMomentumG4;
    TH1F* _hStepLength;
    
    TH1F* _hEdep;
    TH1F* _hEdepMC;
    TH1F* _hNcrystal;
    TH1F* _hNcrstep;
    TH1F* _hNrostep;
    TH1F* _hEdepROMC;

    TH1F* _hRCEdep;
    TH1F* _hRCTime;
    TH1F* _hRCNCrystals;

    TH1F* _hRCEdepMC;
    TH1F* _hRCTimeMC;
    TH1F* _hRCNCrystalsMC;

    TH1F* _hTargetEdep;
    TH1F* _hTargetPathLength;
    TH1F* _hTargetNfoils;
    TH2F* _hTargetNfoils2D;

    TNtuple* _ntup;
    TGraph*  _xyHits;

    // Need to keep track of TGraph entries by hand.
    int _xyHitCount;

    // CRV
    TH1F*    _hCRVMultiplicity;
    TNtuple* _ntupCRV;

    int _nBadG4Status;

    // Do the work specific to one of the trackers.
    void doLTracker(const edm::Event& event);
    void doITracker(const edm::Event& event);
    void doCalorimeter(const edm::Event& event);
    void doStoppingTarget(const edm::Event& event);
    void doCRV(const edm::Event& event);

    // A helper function.
    int countHitNeighbours( Straw const& straw, 
                            edm::Handle<StepPointMCCollection>& hits );

  };
  
} // end namespace mu2e
