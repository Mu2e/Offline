//
// Plugin to test that I can read back the persistent data about straw hits.
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: ReadStrawHit_module.cc,v 1.14 2012/08/22 22:21:32 genser Exp $
// $Author: genser $
// $Date: 2012/08/22 22:21:32 $
//
// Original author Rob Kutschke. Updated by Ivan Logashenko.
//                               Updated by KLG
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Provenance.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"

// CLHEP includes.
#include "CLHEP/Random/RandGaussQ.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "ConditionsService/inc/MassCache.hh"
#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/StepPointMCStrawHit.hh"
#include "HitMakers/inc/formStepPointMCStrawHit.hh"


using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class ReadStrawHit : public art::EDAnalyzer {
  public:
    explicit ReadStrawHit(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _minimumLength(pset.get<double>("minimumLength",0.01)),   // mm
      _enableFlightTimeCorrection(pset.get<bool>("flightTimeCorrection",false)),
      _gaussian( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
      _hHitTime(0),
      _hHitDeltaTime(0),
      _hHitAmplitude(0),
      _hHitEnergy(0),
      _hNHits(0),
      _hNHitsPerWire(0),
      _hDriftTime(0),
      _hDriftDistance(0),
      _hDistanceToMid(0),
      _hNG4Steps(0),
      _hG4StepLength(0),
      _hG4StepEdep(0),
      _ntup(0),
      _detntup(0)
    {
    }
    virtual ~ReadStrawHit() { }

    virtual void beginJob();

    void analyze( art::Event const& e);

  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

    // Some reconstruction related data
    double _minimumLength;
    double _enableFlightTimeCorrection;
    // Random number distributions
    CLHEP::RandGaussQ _gaussian;

    // Some diagnostic histograms.
    TH1F* _hHitTime;
    TH1F* _hHitDeltaTime;
    TH1F* _hHitAmplitude;
    TH1F* _hHitEnergy;
    TH1F* _hNHits;
    TH1F* _hNHitsPerWire;
    TH1F* _hDriftTime;
    TH1F* _hDriftDistance;
    TH1F* _hDistanceToMid;
    TH1F* _hNG4Steps;
    TH1F* _hG4StepLength;
    TH1F* _hG4StepEdep;
    TNtuple* _ntup;
    TNtuple* _detntup;

  };

  void ReadStrawHit::beginJob(){


    if ( _diagLevel > 0 ) {
      cout << "ReadStrawHit Diaglevel: "
           << _diagLevel << " "
           << _maxFullPrint
           << endl;
    }

    art::ServiceHandle<art::TFileService> tfs;

    _hHitTime      = tfs->make<TH1F>( "hHitTime",      "Hit Time (ns)", 200, 0., 2000. );
    _hHitDeltaTime = tfs->make<TH1F>( "hHitDeltaTime", "Hit Delta Time (ns)", 80, -20.0, 20. );
    //    _hHitAmplitude = tfs->make<TH1F>( "hHitAmplitude", "Hit Amplitudes (uV)",  100, 0., 100. );
    _hHitEnergy    = tfs->make<TH1F>( "hHitEnergy",    "Hit Energy (keV)", 100, 0., 100. );
    _hNHits        = tfs->make<TH1F>( "hNHits",        "Number of straw hits", 500, 0., 10000. );
    _hNHitsPerWire = tfs->make<TH1F>( "hNHitsPerWire", "Number of hits per straw", 10, 0., 10. );
    _hDriftTime    = tfs->make<TH1F>( "hDriftTime",    "Drift time, ns", 100, 0., 100. );
    _hDriftDistance= tfs->make<TH1F>( "hDriftDistance","Drift Distance, mm", 100, 0., 3. );
    _hDistanceToMid= tfs->make<TH1F>( "hDistanceToMid","Distance to wire center, mm", 160, -1600., 1600. );
    _hNG4Steps     = tfs->make<TH1F>( "hNG4Steps",     "Number of G4Steps per hit", 100, 0., 100. );
    _hG4StepLength = tfs->make<TH1F>( "hG4StepLength", "Length of G4Steps, mm", 100, 0., 10. );
    _hG4StepEdep   = tfs->make<TH1F>( "hG4StepEdep",   "Energy deposition of G4Steps, keV", 100, 0., 10. );
    _ntup          = tfs->make<TNtuple>( "ntup", "Straw Hit ntuple",
                                         "evt:lay:did:sec:hl:mpx:mpy:mpz:dirx:diry:dirz:time:dtime:eDep:driftT:driftDistance:distanceToMid:id");
    _detntup          = tfs->make<TNtuple>( "detntup", "Straw ntuple",
                                            "id:lay:did:sec:hl:mpx:mpy:mpz:dirx:diry:dirz");
  }

  void
  ReadStrawHit::analyze(art::Event const& evt) {

    static int ncalls(0);
    ++ncalls;

    // Handle to the conditions service
    ConditionsHandle<TrackerCalibrations> trackerCalibrations("ignored");

    /*
    // Print the content of current event
    std::vector<art::Provenance const*> edata;
    evt.getAllProvenance(edata);
    cout << "Event info: "
    << evt.id().event() <<  " "
    << " found " << edata.size() << " objects "
    << endl;
    for( int i=0; i<edata.size(); i++ ) cout << *(edata[i]) << endl;
    */

    // Geometry info for the LTracker.
    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    // Get the persistent data about the StrawHits.

    if (ncalls==1){
      const std::deque<Straw>& allstraws = tracker.getAllStraws();
      float detnt[11];
      for (size_t i = 0;i<allstraws.size();i++)
        {
          Straw str = allstraws[i];
          StrawId sid = str.id();
          LayerId lid = sid.getLayerId();
          DeviceId did = sid.getDeviceId();
          SectorId secid = sid.getSectorId();

          //cout << "index: "  << i << " Layer: "<< lid.getLayer()<< " Device: "<< did <<"  Sector:  "<<secid.getSector()<<endl;
          //cout<<str.getHalfLength()<<endl;
          const CLHEP::Hep3Vector vec3j = str.getMidPoint();
          const CLHEP::Hep3Vector vec3j1 = str.getDirection();
          /*
            cout << i <<
            ","<<lid.getLayer()<<
            ","<<did <<
            ","<<secid.getSector()<<
            ","<<str.getHalfLength()<<
            ","<<vec3j.getX()<<
            ","<<vec3j.getY()<<
            ","<<vec3j.getZ()<<
            ","<<vec3j1.getX()<<
            ","<<vec3j1.getY()<<
            ","<<vec3j1.getZ()<<
            endl;
          */
          // Fill the ntuple.
          detnt[0]  = i;
          detnt[1]  = lid.getLayer();
          detnt[2]  = did;
          detnt[3]  = secid.getSector();
          detnt[4]  = str.getHalfLength();
          detnt[5]  = vec3j.getX();
          detnt[6]  = vec3j.getY();
          detnt[7]  = vec3j.getZ();
          detnt[8]  = vec3j1.getX();
          detnt[9]  = vec3j1.getY();
          detnt[10] = vec3j1.getZ();
          _detntup->Fill(detnt);
        }
    }
    art::Handle<StrawHitCollection> pdataHandle;
    evt.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const& hits = *pdataHandle;

    // Get the persistent data about the StrawHitsMCTruth.
    art::Handle<StrawHitMCTruthCollection> truthHandle;
    evt.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const& hits_truth = *truthHandle;

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const& hits_mcptr = *mcptrHandle;

    // Fill histograms

    _hNHits->Fill(hits.size());

    std::map<StrawIndex,int> nhperwire;

    if ( _diagLevel > 2 ) {
      cout << "ReadStrawHit: Total number of straw hits = " << hits.size() << endl;
    }

    // Cache of recently used masses from the particle data table.
    MassCache cache;

    for ( size_t i=0; i<hits.size(); ++i ) {

      // Access data
      StrawHit             const& hit(hits.at(i));
      StrawHitMCTruth      const& truth(hits_truth.at(i));
      PtrStepPointMCVector const& mcptr(hits_mcptr.at(i));

      // Use data from hits
      _hHitTime->Fill(hit.time());
      _hHitDeltaTime->Fill(hit.dt());
      _hHitEnergy->Fill(hit.energyDep()*1000.0);

      StrawIndex si = hit.strawIndex();

      // Use data from G4 hits
      _hNG4Steps->Fill(mcptr.size());
      for( size_t j=0; j<mcptr.size(); ++j ) {
        StepPointMC const& mchit = *mcptr.at(j);
        _hG4StepLength->Fill(mchit.stepLength());
        if (mchit.strawIndex()!=si) {
          // FIXME: it is an approximation; We plot the "crosstalk
          // edep" which is in principle calculated using amplitudes
          // we may also need to store the crosstalk value in a "truth" object
          _hG4StepEdep->Fill(mchit.eDep()*1000.0*
                             trackerCalibrations->CrossTalk(si,mchit.strawIndex()));
        } else {
          _hG4StepEdep->Fill(mchit.eDep()*1000.0);
        }

      }

      const  unsigned id = si.asUint();
      Straw str = tracker.getStraw(si);
      StrawId sid = str.id();
      LayerId lid = sid.getLayerId();
      DeviceId did = sid.getDeviceId();
      SectorId secid = sid.getSectorId();


      if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {

        //if we were to recalculate some of the quantities using
        //formStepPointMCStrawHit and the first StepPointMC here is how:

        std::auto_ptr<StepPointMCStrawHit> spmcshp = 
          formStepPointMCStrawHit(
                                  mcptr.at(0),
                                  si,
                                  _minimumLength,
                                  _enableFlightTimeCorrection,
                                  cache,
                                  _gaussian,
                                  tracker,
                                  trackerCalibrations);
        cout << "ReadStrawHit: StepPointMCStrawHit # (" << spmcshp->_ptr.id() << 
          " " << spmcshp->_ptr.key() << ")"
             << " DCA=" << spmcshp->_dca
             << " driftTNonSm=" << spmcshp->_driftTimeNonSm
             << " driftT=" << spmcshp->_driftTime
             << " distToMid=" << spmcshp->_distanceToMid
             << " t1=" << spmcshp->_t1
             << " t2=" << spmcshp->_t2
             << " edep=" << spmcshp->_edep
             << endl;

        //       assert (spmcshp->_driftTNonSm == truth.driftTime());
        //       assert (spmcshp->_dca == truth.driftDistance());
        //       assert (spmcshp->_toMid == truth.distanceToMid());
        
      }

      // Use MC truth data
      _hDriftTime->Fill(truth.driftTime());
      _hDriftDistance->Fill(truth.driftDistance());
      _hDistanceToMid->Fill(truth.distanceToMid());

      float nt[18];
      const CLHEP::Hep3Vector smidp  = str.getMidPoint();
      const CLHEP::Hep3Vector sdir   = str.getDirection();
      // Fill the ntuple:
      nt[0]  = evt.id().event();
      nt[1]  = lid.getLayer();
      nt[2]  = did;
      nt[3]  = secid.getSector();
      nt[4]  = str.getHalfLength();
      nt[5]  = smidp.getX();
      nt[6]  = smidp.getY();
      nt[7]  = smidp.getZ();
      nt[8]  = sdir.getX();
      nt[9]  = sdir.getY();
      nt[10] = sdir.getZ();
      nt[11] = hit.time();
      nt[12] = hit.dt();
      nt[13] = hit.energyDep();
      nt[14] = truth.driftTime();
      nt[15] = truth.driftDistance();
      nt[16] = truth.distanceToMid();
      nt[17] = id;
      _ntup->Fill(nt);
      // Calculate number of hits per wire
      ++nhperwire[hit.strawIndex()];

      if ( int(evt.id().event()) < _maxFullPrint ) {
        cout << "ReadStrawHit: " 
             << evt.id().event()      << " #" 
             << si                    << " "
             << sid                   << " "
             << hit.time()            << " " 
             << hit.dt()              << " " 
             << hit.energyDep()       << " "
             << truth.driftTime()     << " "
             << truth.driftDistance() << " |";
        for ( size_t j=0; j<mcptr.size(); ++j ) {
          StepPointMC const& mchit = *mcptr.at(j);
          cout << " " << mchit.ionizingEdep();
        }
        cout << endl;
        if ( _diagLevel > 2 ) {
          for( size_t j=0; j<mcptr.size(); ++j ) {
            StepPointMC const& mchit = *mcptr.at(j);
            cout << "ReadStrawHit: StepPointMC #" << j << " : length=" << mchit.stepLength()
                 << " energy=" << mchit.totalEDep() << " time=" <<  mchit.time()
                 << endl;
            // we'll print the track info

            SimParticle sp = *mchit.simParticle();

            cout << "ReadStrawHit: mchit.simParticle() id, pdgid, parent id, pdgid : " 
                 << sp.id() << ", " << sp.pdgId();
            
            bool isSec = sp.isSecondary();

            if ( isSec ) {

              while ( isSec ) {

                SimParticle parent = *sp.parent();

                cout << " is a secondary of "  << parent.id() 
                     << ", " << parent.pdgId();
              
                isSec = parent.isSecondary();
                sp = parent;

              }

              cout << " which is a primary" << endl;

            }
            else {

              cout << " is a primary" << endl;

            }

          }

        }

      }

    }

    for( std::map<StrawIndex,int>::iterator it=nhperwire.begin(); it!= nhperwire.end(); ++it ) {
      _hNHitsPerWire->Fill(it->second);
    }

  } // end of ::analyze.

}


using mu2e::ReadStrawHit;
DEFINE_ART_MODULE(ReadStrawHit);
