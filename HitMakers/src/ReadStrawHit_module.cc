//
// Plugin to test that I can read back the persistent data about straw hits.
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
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
#include "canvas/Utilities/Exception.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class ReadStrawHit : public art::EDAnalyzer {
  public:
    explicit ReadStrawHit(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;

    void analyze( art::Event const& e) override;

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

    // Some diagnostic histograms.
    TH1F* _hHitTime;
    TH1F* _hHitTime1;
    TH1F* _hHitDeltaTime;
    TH1F* _hHitAmplitude;
    TH1F* _hHitEnergy;
    TH1F* _hNHits;
    TH1F* _hNHits1;
    TH1F* _hNHitsPerWire;
    TH1F* _hDriftTime;
    TH1F* _hDriftDistance;
    TH1F* _hDistanceToMid;
    TH1F* _hFractionalDistanceToMid;
    TH1F* _hNG4Steps;
    TH1F* _hG4StepLength;
    TH1F* _hG4StepEdep;
    TH1F* _hG4StepRelTimes;
    TNtuple* _ntup;
    TNtuple* _detntup;

  };

  ReadStrawHit::ReadStrawHit(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _maxFullPrint(pset.get<int>("maxFullPrint",5)),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
    _minimumLength(pset.get<double>("minimumLength",0.01)),   // mm
    _enableFlightTimeCorrection(pset.get<bool>("flightTimeCorrection",false)),
    _hHitTime(0),
    _hHitTime1(0),
    _hHitDeltaTime(0),
    _hHitAmplitude(0),
    _hHitEnergy(0),
    _hNHits(0),
    _hNHits1(0),
    _hNHitsPerWire(0),
    _hDriftTime(0),
    _hDriftDistance(0),
    _hDistanceToMid(0),
    _hFractionalDistanceToMid(0),
    _hNG4Steps(0),
    _hG4StepLength(0),
    _hG4StepEdep(0),
    _hG4StepRelTimes(0),
    _ntup(0),
    _detntup(0){
  }

  void ReadStrawHit::beginJob(){

    if ( _diagLevel > 0 ) {
      cout << "ReadStrawHit Diaglevel: "
           << _diagLevel << " "
           << _maxFullPrint
           << endl;
    }

    art::ServiceHandle<art::TFileService> tfs;

    _hHitTime      = tfs->make<TH1F>( "hHitTime",      "Hit Time (ns)", 200, 0.,  2000. );
    _hHitTime1     = tfs->make<TH1F>( "hHitTime1",     "Hit Time (ns)", 200, 0., 20000. );
    _hHitDeltaTime = tfs->make<TH1F>( "hHitDeltaTime", "Hit Delta Time (ns)", 80, -20.0, 20. );
    //    _hHitAmplitude = tfs->make<TH1F>( "hHitAmplitude", "Hit Amplitudes (uV)",  100, 0., 100. );
    _hHitEnergy    = tfs->make<TH1F>( "hHitEnergy",    "Hit Energy (keV)", 100, 0., 100. );
    _hNHits        = tfs->make<TH1F>( "hNHits",        "Number of straw hits", 500, 0.,  500. );
    _hNHits1       = tfs->make<TH1F>( "hNHits1",       "Number of straw hits", 200, 0., 4000. );
    _hNHitsPerWire = tfs->make<TH1F>( "hNHitsPerWire", "Number of hits per straw", 10, 0., 10. );
    _hDriftTime    = tfs->make<TH1F>( "hDriftTime",    "Drift time, ns", 100, 0., 100. );
    _hDriftDistance= tfs->make<TH1F>( "hDriftDistance","Drift Distance, mm", 100, 0., 3. );
    _hDistanceToMid= tfs->make<TH1F>( "hDistanceToMid","Distance to wire center, mm", 160, -1600., 1600. );
    _hFractionalDistanceToMid
                   = tfs->make<TH1F>( "hFractionalDistanceToMid","(Distance to wire center)(half length of wire)", 100, -1., 1. );
    _hNG4Steps     = tfs->make<TH1F>( "hNG4Steps",     "Number of G4Steps per hit", 100, 0., 100. );
    _hG4StepLength = tfs->make<TH1F>( "hG4StepLength", "Length of G4Steps, mm", 100, 0., 10. );
    _hG4StepEdep   = tfs->make<TH1F>( "hG4StepEdep",   "Energy deposition of G4Steps, keV", 100, 0., 10. );
    _hG4StepRelTimes = tfs->make<TH1F>( "hG4StepRelTimes", "Hit Relative Times of G4Steps, ns", 100, -100., 100.);

    _ntup          = tfs->make<TNtuple>( "ntup", "Straw Hit ntuple",
                                         "evt:lay:did:sec:hl:mpx:mpy:mpz:dirx:diry:dirz:time:dtime:eDep:driftT:driftDistance:distanceToMid:id:hitx:hity:hitz:fracDistanceToMid");
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
    cout << __func__ << " Event info: "
    << evt.id().event() <<  " "
    << " found " << edata.size() << " objects "
    << endl;
    for( size_t i=0; i<edata.size(); i++ ) cout << *(edata[i]) << endl;
    */

    // Geometry info for the Tracker.
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
          PlaneId did = sid.getPlaneId();
          PanelId secid = sid.getPanelId();

          // cout <<  __func__ << " index: "  << i << " Layer: "<< lid.getLayer()<< " Plane: "<< did <<"  Panel:  "<<secid.getPanel()<<endl;
          // cout<<str.getHalfLength()<<endl;
          const CLHEP::Hep3Vector vec3j = str.getMidPoint();
          const CLHEP::Hep3Vector vec3j1 = str.getDirection();
          /*
            cout << i <<
            ","<<lid.getLayer()<<
            ","<<did <<
            ","<<secid.getPanel()<<
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
          detnt[3]  = secid.getPanel();
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
    bool gblresult;

    art::Handle<StrawHitCollection> pdataHandle;
    gblresult = evt.getByLabel(_makerModuleLabel,pdataHandle);
    ( _diagLevel > 0 ) &&
      std::cout
      << __func__ << " getting data by getByLabel: label, instance, result " << std::endl
      << " StrawHitCollection             _makerModuleLabel: " << _makerModuleLabel << ", "
      << ", " << gblresult << std::endl;

    if (!gblresult) throw cet::exception("DATA") << " Missing data";

    StrawHitCollection const& hits = *pdataHandle;

    // Get the persistent data about the StrawHitsMCTruth or StrawDigiMC
    // fixme: we may templetize this
    art::Handle<StrawHitMCTruthCollection> truthHitMCHandle;
    bool gblSHresult = evt.getByLabel(_makerModuleLabel,truthHitMCHandle);
    ( _diagLevel > 0 ) &&
      std::cout
      << __func__ << " getting data by getByLabel: label, instance, result " << std::endl
      << " StrawHitMCTruthCollection      _makerModuleLabel: " << _makerModuleLabel << ", "
      << ", " << gblresult << std::endl;

    StrawHitMCTruthCollection const& hits_SHtruth = gblSHresult ?
      *truthHitMCHandle : StrawHitMCTruthCollection();

    art::Handle<StrawDigiMCCollection>     truthDigiMCHandle;
    bool gblSDresult = evt.getByLabel(_makerModuleLabel,truthDigiMCHandle);
    ( _diagLevel > 0 ) && 
      std::cout 
      << __func__ << " getting data by getByLabel: label, instance, result " << std::endl
      << " StrawHitMCTruthCollection      _makerModuleLabel: " << _makerModuleLabel << ", StrawHitMC"
      << ", " << gblresult << std::endl;

    StrawDigiMCCollection const& hits_SDtruth = gblSDresult ?
      *truthDigiMCHandle : StrawDigiMCCollection();

    if (!(gblSHresult||gblSDresult)) throw cet::exception("DATA") << " Missing data";

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    gblresult = evt.getByLabel(_makerModuleLabel,mcptrHandle);

    ( _diagLevel > 0 ) &&
      std::cout
      << __func__ << " getting data by getByLabel: label, instance, result " << std::endl
      << " PtrStepPointMCVectorCollection _makerModuleLabel: " << _makerModuleLabel << ", StrawHitMCPtr"
      << ", " << gblresult << std::endl;

    if (!gblresult) throw cet::exception("DATA") << " Missing data";

    PtrStepPointMCVectorCollection const& hits_mcptr = *mcptrHandle;

    // Fill histograms

    _hNHits->Fill(hits.size());
    _hNHits1->Fill(hits.size());

    std::map<StrawIndex,int> nhperwire;

    if ( _diagLevel > 2 ) {
      cout << "ReadStrawHit: Total number of straw hits = " << hits.size() << endl;
    }

    for ( size_t i=0; i<hits.size(); ++i ) {

      // Access data
      StrawHit             const& hit(hits.at(i));
      StrawHitMCTruth      const& SHtruth = gblSHresult ? hits_SHtruth.at(i) : StrawHitMCTruth() ;
      StrawDigiMC          const& SDtruth = gblSDresult ? hits_SDtruth.at(i) : StrawDigiMC() ;
      PtrStepPointMCVector const& mcptr(hits_mcptr.at(i));

      // Use data from hits
      _hHitTime ->Fill(hit.time());
      _hHitTime1->Fill(hit.time());
      _hHitDeltaTime->Fill(hit.dt());
      _hHitEnergy->Fill(hit.energyDep()*1000.0);

      StrawIndex si = hit.strawIndex();

      // Use data from G4 hits
      _hNG4Steps->Fill(mcptr.size());
      for( size_t j=0; j<mcptr.size(); ++j ) {
        StepPointMC const& mchit = *mcptr.at(j);
        _hG4StepLength->Fill(mchit.stepLength());
        //        _hG4StepRelTimes->Fill(fabs(mchit.time()-hit.time()));
        _hG4StepRelTimes->Fill((mchit.time()-hit.time()));
        // step time rel to the formed hit time; FIXME fabs should not be needed
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
      PlaneId did = sid.getPlaneId();
      PanelId secid = sid.getPanelId();

      double fracDist = 0.0;
      StrawEnd itdc(TrkTypes::cal);

      if (gblSHresult) {
        fracDist = SHtruth.distanceToMid()/str.getHalfLength();
      } else {
	fracDist = SDtruth.distanceToMid(itdc)/str.getHalfLength();
      }

      // Use MC truth data

      if (gblSHresult) {
        _hDriftTime->Fill(SHtruth.driftTime());
        _hDriftDistance->Fill(SHtruth.driftDistance());
        _hDistanceToMid->Fill(SHtruth.distanceToMid());
        _hFractionalDistanceToMid->Fill(fracDist);
      } else {
        _hDriftTime->Fill(0.0);
        _hDriftDistance->Fill(SDtruth.driftDistance(itdc));
        _hDistanceToMid->Fill(SDtruth.distanceToMid(itdc));
        _hFractionalDistanceToMid->Fill(fracDist);
      }

      const CLHEP::Hep3Vector smidp  = str.getMidPoint();
      const CLHEP::Hep3Vector sdir   = str.getDirection();

      // calculate the hit position
      SHInfo strawHitInfo;
      trackerCalibrations->StrawHitInfo(str, hit, strawHitInfo);

      // we may also need the truth hit position, do we need another function?

      // Fill the ntuple:
      float nt[_ntup->GetNvar()];
      nt[0]  = evt.id().event();
      nt[1]  = lid.getLayer();
      nt[2]  = did;
      nt[3]  = secid.getPanel();
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
      nt[14] = gblSHresult ? SHtruth.driftTime() : 0.0;
      nt[15] = gblSHresult ? SHtruth.driftDistance() : SDtruth.driftDistance(itdc);
      nt[16] = gblSHresult ? SHtruth.distanceToMid() : SDtruth.distanceToMid(itdc);
      nt[17] = id;
      nt[18] = strawHitInfo._pos.getX();
      nt[19] = strawHitInfo._pos.getY();
      nt[20] = strawHitInfo._pos.getZ();
      nt[21] = fracDist;

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
             << ( gblSHresult ? SHtruth.driftTime() : 0.0 ) << " "
             << ( gblSHresult ? SHtruth.driftDistance() : SDtruth.distanceToMid(itdc) ) << " |";
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
