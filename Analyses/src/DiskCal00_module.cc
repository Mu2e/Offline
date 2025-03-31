//
// Check self consistency of hits in the Disk Calorimeter.
//
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>
#include <set>
#include <cmath>
#include <algorithm>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

// Mu2e includes.
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"

// Root includes.
#include "TH1F.h"

// Other includes.
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  class DiskCal00 : public art::EDAnalyzer {
  public:

    explicit DiskCal00(fhicl::ParameterSet const& pset);

    virtual void beginJob( );
    virtual void endJob  ( );

    virtual void beginRun( const art::Run& run);
    virtual void endRun  ( const art::Run& run);

    virtual void analyze ( const art::Event& event);


  private:

    // Module label of the g4 module that made the hits.
    std::string _apdHitMakerModuleLabel;
    std::string _crystalHitMakerModuleLabel;

    const PtrStepPointMCVectorCollection*   steps;

    int    maxId;
    size_t maxHits;
    size_t maxCrystalsHit;

    TH1F* _hNHits;
    TH1F* _hNCrystal;
    TH1F* _hTotalEDep;
    TH1F* _hEnergyDep;
    TH1F* _hEnergyDep1;
    TH1F* _hEnergyDep2;
    TH1F* _hEnergyDep3;
    TH1F* _hTime;
    TH1F* _hStepTime;
    TH1F* _hNSiPMID;
    TH1F* _hNSteps;
    TH1F* _hDeltaTime;

    void printCalInfo();
    void followHistory( const CaloHit& cryHit);
    void simCheck(const art::Event& event);

  };

  DiskCal00::DiskCal00(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    _apdHitMakerModuleLabel("CaloReadoutHitsMaker"),
    _crystalHitMakerModuleLabel("CaloHitsMaker"),

    //
    steps(0),

    // Limit checks
    maxId(0),
    maxHits(0),
    maxCrystalsHit(0),

    // Histograms
    _hNHits(0),
    _hNCrystal(0),
    _hTotalEDep(0),
    _hEnergyDep(0),
    _hEnergyDep1(0),
    _hEnergyDep2(0),
    _hEnergyDep3(0),
    _hTime(0),
    _hStepTime(0),
    _hNSiPMID(0),
    _hNSteps(0),
    _hDeltaTime(0)
  {}

  void DiskCal00::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _hNHits       = tfs->make<TH1F>( "hNHits",      "Number of hit in the event",                     100,  0.,  4000. );
    _hNCrystal    = tfs->make<TH1F>( "hNcrystal",   "Number of hit crystals in calorimeter",          100,  0.,  2000. );
    _hTotalEDep   = tfs->make<TH1F>( "hTotalEDep",  "Total energy deposition in calorimeter;[MeV]",   100,  0.,  3000. );
    _hEnergyDep   = tfs->make<TH1F>( "hEnergyDep",  "Energy deposit in one Crystal;[MeV}",            100,  0.,   100. );
    _hEnergyDep1  = tfs->make<TH1F>( "hEnergyDep1", "Energy deposit in one Crystal;[MeV]",            100,  0.,     1. );
    _hEnergyDep2  = tfs->make<TH1F>( "hEnergyDep2", "Log10(Energy deposit in one Crystal, MeV)",      180, -6.,     3. );
    _hEnergyDep3  = tfs->make<TH1F>( "hEnergyDep3", "Log10(TotalEnergy deposit in one Crystal, Mev)", 180, -6.,     3. );
    _hTime        = tfs->make<TH1F>( "hTime",       "Time of hit;[ns]",                               125,  0.,  2500. );
    _hStepTime    = tfs->make<TH1F>( "hStepTime",   "Time of hit;[ns]",                               125,  0.,  2500. );
    _hNSiPMID     = tfs->make<TH1F>( "hNSiPMID",    "Number of SiPM Ids",                               4, -1.,     3. );
    _hNSteps      = tfs->make<TH1F>( "hNSteps",     "Number of StepPointMCs per Crystal",             100,  0.,   500. );
    _hDeltaTime   = tfs->make<TH1F>( "hDeltaTime",  "Delta Time (step-hit)",                          110, -10.,  550. );

  } // end beginJob

  void DiskCal00::endJob(){
    cout << "Largest Id encountered:          " << maxId   << endl;
    cout << "Maximum number of  hits:         " << maxHits << endl;
    cout << "Maximum number of hits crystals: " << maxCrystalsHit << endl;
  }

  void DiskCal00::beginRun( const art::Run& ){
  }

  void DiskCal00::endRun( const art::Run& ){
    printCalInfo();
  }

  void DiskCal00::analyze(const art::Event& event) {

    simCheck(event);

    art::Handle<CaloHitCollection> cryHitsHandle;
    event.getByLabel(_crystalHitMakerModuleLabel,cryHitsHandle);
    CaloHitCollection const& cryHits = *cryHitsHandle;

    art::Handle<PtrStepPointMCVectorCollection> stepsHandle;
    event.getByLabel(_apdHitMakerModuleLabel, "CaloHitMCCrystalPtr", stepsHandle);
    steps = stepsHandle.product();

    _hNHits->Fill(cryHits.size());
    maxHits = std::max( cryHits.size(), maxHits );
    GeomHandle<DiskCalorimeter> calGeom;

    double totalEdep(0.);
    set<int> hit_crystals;
    set<int> crystalsOverPed;

    double minEnergy(0.001);

    for ( CaloHitCollection::const_iterator i=cryHits.begin(), e=cryHits.end(); i != e; ++i ){

      const CaloHit& cryHit(*i);

      totalEdep += cryHit.energyDep();
      _hEnergyDep ->Fill(cryHit.energyDep());
      _hEnergyDep1->Fill(cryHit.energyDep());
      double logE    = ( cryHit.energyDep() > 0.)      ? log10(cryHit.energyDep())      : -5.5;
      double logETot = ( cryHit.energyDepTot() > 0.) ? log10(cryHit.energyDepTot()) : -5.5;
      _hEnergyDep2->Fill(logE);
      _hEnergyDep3->Fill(logETot);
      _hTime->Fill(cryHit.time());
      _hNSiPMID->Fill( cryHit.nSiPMs() );

      int id = cryHit.crystalID();
      hit_crystals.insert(id);

      maxId = std::max(id, maxId);

      if ( cryHit.energyDep() > minEnergy ){
        crystalsOverPed.insert(id);
      }

      followHistory( cryHit );
    }

    _hTotalEDep->Fill(totalEdep);
    _hNCrystal->Fill(hit_crystals.size());

    maxCrystalsHit = std::max( hit_crystals.size(), maxCrystalsHit );

    int highestId = ( hit_crystals.empty() ) ? -1 : *hit_crystals.rbegin();

    static int nPrint(0);
    if ( ++nPrint<25 ){
      cout << "Sizes: "
           << cryHits.size()         << " "
           << hit_crystals.size() << " "
           << crystalsOverPed.size() << " "
           << totalEdep << " "
           << highestId << " "
           << steps->size()
           << endl;
    }

  } // end analyze

  void DiskCal00::followHistory( const CaloHit& cryHit){

    /*
    vector<art::Ptr<CaloHit> > const & readouts(cryHit.readouts());

    if ( readouts.size()< 2 ) return;

    for ( vector<art::Ptr<CaloHit> >::const_iterator i=readouts.begin(), e=readouts.end(); i != e; ++i ){
      vector<art::Ptr<StepPointMC> > const& stepv = steps->at(i->key());

      _hNSteps->Fill(stepv.size());

      for (  vector<art::Ptr<StepPointMC> >::const_iterator j=stepv.begin(), je=stepv.end(); j != je; ++j ){
        const StepPointMC& step = **j;
        if ( cryHit.crystalID() != int(step.volumeId()) ){
          cout << "FUBAR: " << cryHit.crystalID() << " " << step.volumeId() << endl;
        }
        _hStepTime->Fill(step.time());
        _hDeltaTime->Fill(step.time()-cryHit.time());
      }
    }
    */
  }

  void DiskCal00::simCheck(const art::Event& event) {
    typedef std::vector< art::Handle<SimParticleCollection> > HandleVector;
    HandleVector allSims = event.getMany<SimParticleCollection>();
    cout << "Size: " << allSims.size() << endl;

    for ( HandleVector::const_iterator i=allSims.begin(), e=allSims.end(); i != e; ++i ){
      const SimParticleCollection& sims(**i);
      for ( SimParticleCollection::const_iterator j=sims.begin(), je=sims.end(); j != je; ++j ){
        const SimParticle& sim = j->second;
        if ( sim.isSecondary() ){
          double t  = sim.startGlobalTime();
          double t0 =  sim.parent()->startGlobalTime();
          //double dt = t-t0;
          if ( t < t0 ){
            cout << "FUBAR: " << i-allSims.begin() << " "
                 << t << " " << t0 << " " << sim.id() << " " << sim.parent()->id() << " "
                 << sim.creationCode()
                 << endl;
          } else{
            cout << "OK:    " << i-allSims.begin() << " "
                 << t << " " << t0 << " " << sim.id() << " " << sim.parent()->id() << " "
                 << sim.creationCode()
                 << endl;
          }
        }
      }
    }

  }

  void DiskCal00::printCalInfo(){
    DiskCalorimeter const& cal(*GeomHandle<DiskCalorimeter>());
    int nSiPM = cal.nCrystals()*cal.caloInfo().getInt("nSiPMPerCrystal");
    cout << "Information about the disk Calorimeter: "  << endl;
    cout << "Number of disks:    " << cal.nDisks()      << endl;
    cout << "Number of Readouts: " << nSiPM << " "  << CaloConst::_nSiPMPerCrystal << " " << nSiPM/CaloConst::_nSiPMPerCrystal << endl;
    cout << "Hex side size:      " << 2.0*cal.caloInfo().getDouble("crystalXYLength") << endl;

    cout << "Depth:              " << cal.caloInfo().getDouble("crystalZLength")   << endl;
    cout << "Origin:             " << cal.geomUtil().origin()      << endl;
    for (unsigned i=0; i<cal.nDisks(); ++i){
      Disk const& disk = cal.disk(i);
      cout << "Disk: " << i << " " << "origin: " << disk.geomInfo().origin() << endl;
    }
  }

}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
DEFINE_ART_MODULE(mu2e::DiskCal00)
