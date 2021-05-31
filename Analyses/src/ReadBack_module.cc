//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// Original author Rob Kutschke
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBarDetail.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  class ReadBack : public art::EDAnalyzer {
  public:

    explicit ReadBack(fhicl::ParameterSet const& pset);
    virtual ~ReadBack() { }

    virtual void beginJob() override;
    virtual void endJob() override;

    // This is called for each event.
    virtual void analyze(const art::Event& e) override;

  private:

    typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;

    // Start: run time parameters

    // Diagnostics printout level
    int _diagLevel;

    // Input tag of the generated particles.
    art::InputTag _gensTag;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Module label of the generator module that was passed as input to G4.
    std::string _generatorModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Name of the calorimeter StepPoint collection
    std::string _calorimeterStepPoints;

    // Module which made the MC CaloShowers
    std::string _caloShowerSimModuleLabel;

    // Module which made the CaloHits
    std::string _caloCrystalModuleLabel;

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

    // by how much the straw gas and hit positions can differ
    double _strawHitPositionTolerance;

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

    TH1F* _hCaStepEdep;
    TH1F* _hCaStepNum;
    TH1F* _hCaShowerEdep;
    TH1F* _hCaShowerNum;
    TH1F* _hCaTime;
    TH1F* _hCaEdep;
    TH2F* _hCaCrystalXY;
    TH1F* _hCaEdepMC;
    TH1F* _hCaNcrystal;

    TH1F* _hTargetEdep;
    TH1F* _hTargetPathLength;
    TH1F* _hTargetNfoils;
    TH2F* _hTargetNfoils2D;

    TH1F* _hNPointTrajectories;
    TH1F* _hNPointsPerTrajectory;
    TH1F* _hNPointsPerTrajectoryLog;

    TNtuple* _ntup;
    TGraph*  _xyHits;

    // Need to keep track of TGraph entries by hand.
    int _xyHitCount;

    // CRV
    TH1F*    _hCRVMultiplicity;
    TNtuple* _ntupCRV;

    int _nBadG4Status;

    // Examine various parts of the event.
    void doTracker         ( const art::Event& event );
    void doCalorimeter      ( const art::Event& event );
    void doStoppingTarget   ( const art::Event& event );
    void doCRV              ( const art::Event& event );
    void doPointTrajectories( const art::Event& event );

    // A helper function.
    int countHitNeighbours( Straw const& straw,
                            art::Handle<StepPointMCCollection>& hits );

  };
  ReadBack::ReadBack(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    // Run time parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _gensTag(pset.get<string>("gensTag","generate")),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _generatorModuleLabel(pset.get<string>("generatorModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _calorimeterStepPoints(pset.get<string>("calorimeterStepPoints","calorimeter")),
    _caloShowerSimModuleLabel(pset.get<string>("caloShowerSimModuleLabel","CaloShowerROMaker")),
    _caloCrystalModuleLabel(pset.get<string>("caloCrystalModuleLabel","CaloHitMaker")),
    _targetStepPoints(pset.get<string>("targetStepPoints","stoppingtarget")),
    _crvStepPoints(pset.get<string>("CRVStepPoints","CRV")),
    _minimumEnergy(pset.get<double>("minimumEnergy")),
    _maxFullPrint(pset.get<int>("maxFullPrint",5)),
    _xyHitsMax(pset.get<int>("xyHitsMax",10000)),
    _strawHitPositionTolerance(pset.get<double>("SHPositionTolerance",0.01)),
    // looking for gross errors only

    // Histograms
    _nAnalyzed(0),
    _hRadius(0),
    _hTime(0),
    _hMultiplicity(0),
    _hDriftDist(0),
    _hDriftDistW(0),
    _hxHit(0),
    _hyHit(0),
    _hzHit(0),
    _hHitNeighbours(0),
    _hCheckPointRadius(0),
    _hCheckPointRadiusW(0),
    _hCheckPointWireZ(0),
    _hMomentumG4(0),
    _hStepLength(0),
    _hCaStepEdep(0),
    _hCaStepNum(0),
    _hCaShowerEdep(0),
    _hCaShowerNum(0),
    _hCaTime(0),
    _hCaEdep(0),
    _hCaCrystalXY(0),
    _hCaEdepMC(0),
    _hCaNcrystal(0),
    _hTargetEdep(0),
    _hTargetPathLength(0),
    _hTargetNfoils(0),
    _hTargetNfoils2D(0),
    _hNPointTrajectories(0),
    _hNPointsPerTrajectory(0),
    _hNPointsPerTrajectoryLog(0),
    _ntup(0),
    _xyHits(0),
    _xyHitCount(0),
    // CRV
    _hCRVMultiplicity(0),
    _ntupCRV(0),
    // Remaining member data
    _nBadG4Status(0){
  }

  void ReadBack::beginJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    // Create some 1D histograms.
    _hRadius       = tfs->make<TH1F>( "hRadius",  "Radius of StepPoint Locations;(mm)",     500,  300.,  800. );
    _hEnergyDep    = tfs->make<TH1F>( "hEnergyDep",    "Energy Deposited;(keV)",  100,  0.,   10. );
    _hTime         = tfs->make<TH1F>( "hTime",         "Pulse Height;(ns)",       100,  0., 2000. );
    _hMultiplicity = tfs->make<TH1F>( "hMultiplicity", "StepPoints per Event",         300,  0.,  300. );
    _hDriftDist    = tfs->make<TH1F>( "hDriftDist", "Crude Drift Distance;(mm)",  100,  0.,   3.  );
    _hDriftDistW   = tfs->make<TH1F>( "hDriftDistW", "Normalized Crude Drift Distance", 150,  0., 1.5 );

    _hxHit         = tfs->make<TH1F>( "hxHit",  "X of StepPoint;(mm)",                  1600, -800., 800. );
    _hyHit         = tfs->make<TH1F>( "hyHit",  "Y of StepPoint;(mm)",                  1600, -800., 800. );
    _hzHit         = tfs->make<TH1F>( "hzHit",  "Z of StepPoint;(mm)",                  3000, -1500., 1500. );

    _hHitNeighbours    = tfs->make<TH1F>( "hHitNeighbours",  "Number of hit neighbours",
                                          10, 0., 10. );

    _hCheckPointRadius  = tfs->make<TH1F>( "hCheckPointRadius", "Radius of Reference point; (mm)",
                                           100, 2.25, 2.75 );
    _hCheckPointRadiusW = tfs->make<TH1F>( "hCheckPointRadiusW","Normalized Radius of Reference point",
                                           200, 0., 2. );
    _hCheckPointWireZ   = tfs->make<TH1F>( "hCheckPointWireZ",  "Normalized Wire Z of Reference point",
                                           300, -1.5, 1.5 );

    _hMomentumG4 = tfs->make<TH1F>( "hMomentumG4",  "Mommenta of particles created inside G4; (MeV)",
                                    100, 0., 100. );
    _hStepLength = tfs->make<TH1F>( "hStepLength",  "G4 Step Length in Sensitive Detector; (mm)",
                                    100, 0., 10. );

    _hCaStepEdep =   tfs->make<TH1F>("hCaStepEdep",  "Total energy deposition from StepPt in each crystal",
                                     2000, 0., 2000. );
    _hCaStepNum    = tfs->make<TH1F>("hCaStepNum",   "Total number of StepPt in each crystal",
                                     2000, 0., 2000. );
    _hCaCrystalXY  = tfs->make<TH2F>("hCaCrystalXY", "calorimeter stepPt position in crystal frame",
                                     100,-50.,50.,  100,-50.,50. );

    _hCaShowerEdep = tfs->make<TH1F>("hCaShowerEdep",  "Total energy deposition from StepPt in each crystal",
                                     2000, 0., 2000. );
    _hCaShowerNum  = tfs->make<TH1F>("hCaShowerNum",   "Total number of StepPt in each crystal",
                                     2000, 0., 2000. );
    _hCaTime       = tfs->make<TH1F>("hCaTime",   "Calorimeter hit time",
                                     2000, 0., 2000. );
    _hCaEdep       = tfs->make<TH1F>("hCaEdep",      "Total energy deposition in calorimeter",
                                     200, 0., 200. );
    _hCaNcrystal   = tfs->make<TH1F>("hCaNcrystal",    "Number of crystals hit",
                                     2000,  0., 2000. );

    // Stopping target histograms

    _hTargetEdep = tfs->make<TH1F>( "hTargetEdep",
                                    "Energy deposition in the stopping target",
                                    100, 0., 5. );
    _hTargetPathLength = tfs->make<TH1F>( "hTargetPathLength",
                                          "Path length in the stopping target",
                                          100, 0., 5. );
    _hTargetNfoils = tfs->make<TH1F>( "hTargetNfoils",
                                      "Number of stopping target foils crossed by particle",
                                      20, 0., 20. );
    _hTargetNfoils2D = tfs->make<TH2F>( "hTargetNfoils2D",
                                        "Number of stopping target foils vs foil of origin",
                                        20, 0., 20., 20, 0, 20. );

    // Histograms for the MCTrajectory Collection
    _hNPointTrajectories      = tfs->make<TH1F>( "hNPointTrajectories",      "Number of PointTrajectories stored",      50,  0.,  50. );
    _hNPointsPerTrajectory    = tfs->make<TH1F>( "hNPointsPerTrajectory",    "Number of Points Per Trajectory",        100,  0., 200. );
    _hNPointsPerTrajectoryLog = tfs->make<TH1F>( "hNPointsPerTrajectoryLog", "Log10(Number of Points Per Trajectory)", 120, -1.,   5. );

    // Create tracker ntuple.
    _ntup = tfs->make<TNtuple>( "ntup", "Hit ntuple",
                                "evt:trk:sid:hx:hy:hz:wx:wy:wz:dca:time:pln:pnl:lay:pdgId:genId:edep:p:step:hwz:straw");

    // Create a TGraph;
    // - Syntax to set name and title is weird; that's just root.
    // - Must append to the output file by hand.
    _xyHits = tfs->make<TGraph>(_xyHitsMax);
    _xyHits->SetName("xyHits");
    _xyHits->SetTitle("Y vs X for StepPointMC");
    gDirectory->Append(_xyHits);

    // CRV

    _hCRVMultiplicity = tfs->make<TH1F>( "hCRVMultiplicity", "CRV StepPointMC per Bar", 5000,  0.,  5000. );

    // Create CRV ntuple.
    _ntupCRV = tfs->make<TNtuple>( "ntupCRV", "CRV Hit ntuple",
                                   "evt:trk:bid:hx:hy:hz:bx:by:bz:dx:dy:dz:lx:ly:lz:time:shld:mod:lay:pdgId:genId:edep:p:step");

  }

  void ReadBack::analyze(const art::Event& event) {

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // Inquire about the completion status of G4.
    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( _g4ModuleLabel, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;
    if ( _nAnalyzed < _maxFullPrint ){
      cerr << g4Status << endl;
    }

    // Abort if G4 did not complete correctly.
    // Use your own judgement about whether to abort or to continue.
    if ( g4Status.status() > 1 ) {
      ++_nBadG4Status;
      mf::LogError("G4")
        << "Aborting ReadBack::analyze due to G4 status\n"
        << g4Status;
      return;
    }

    // Call code appropriate for the tracker that is installed in this job.
    art::ServiceHandle<GeometryService> geom;
    if( geom->hasElement<Tracker>() ){
      doTracker(event);
    }

    if (geom->hasElement<Calorimeter>() )doCalorimeter(event);

    doStoppingTarget(event);

    if( geom->hasElement<CosmicRayShield>() )
    {
      doCRV(event);
    }

    doPointTrajectories( event );

  }

  void ReadBack::doCalorimeter(const art::Event& event)
  {

     //Handle to the calorimeter
     art::ServiceHandle<GeometryService> geom;
     Calorimeter const & cal = *(GeomHandle<Calorimeter>());


     // These selectors will select data products with the given instance name, and ignore all other fields of the product ID.
     art::ProductInstanceNameSelector getCrystalSteps(_calorimeterStepPoints);
     HandleVector crystalStepsHandles = event.getMany<StepPointMCCollection>(getCrystalSteps);

     //Calorimeter shower MC
     art::Handle<CaloShowerSimCollection> caloShowerSimHandle;
     event.getByLabel(_caloShowerSimModuleLabel, caloShowerSimHandle);
     const CaloShowerSimCollection& caloShowerSims(*caloShowerSimHandle);

     //Crystal hits (average from readouts)
     art::Handle<CaloHitCollection> CaloHitsHandle;
     event.getByLabel(_caloCrystalModuleLabel, CaloHitsHandle);
     CaloHitCollection const& CaloHits(*CaloHitsHandle);




     //collect all StepPointMC in crystal

     map<int,double> stepPtMap;
     map<int,int>    stepPtMap2;
     for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end(); i != e; ++i )
     {
          const art::Handle<StepPointMCCollection>& handle(*i);
          const StepPointMCCollection& steps(*handle);

          for (const auto& step : steps )
	  {
	      stepPtMap[step.volumeId()] += step.totalEDep();
	      ++stepPtMap2[step.volumeId()];


              CLHEP::Hep3Vector hitPos  = cal.geomUtil().mu2eToCrystal(step.volumeId(),step.position());
	      _hCaCrystalXY->Fill(hitPos.x(),hitPos.y());

              if ( _diagLevel > 1 && _nAnalyzed < _maxFullPrint )
        	cout << "Readback: Calo StepPointMC (" << step.volumeId() << "," << step.totalEDep() << ")";
	  }
     }

     for (const auto& iter : stepPtMap) _hCaStepEdep->Fill(iter.first,iter.second);
     for (const auto& iter : stepPtMap2) _hCaStepNum->Fill(iter.first,iter.second);



     //do the same with the CaloShowers
     if (!caloShowerSimHandle.isValid()) return;

     map<int,double> showerMap;
     map<int,int> showerMap2;
     for (const auto& showerSim : caloShowerSims)
     {
         showerMap[showerSim.crystalID()] += showerSim.energyDep();
         for (const auto& step : showerSim.caloShowerSteps()) showerMap2[showerSim.crystalID()] += step->nCompress();
         _hCaTime->Fill(showerSim.time());

         if ( _diagLevel > 1 && _nAnalyzed < _maxFullPrint )
	   std::cout<<"Readback: caloshower in crystal "<< showerSim.crystalID()<<" eDep = "<<showerSim.energyDep()
	            <<" time = "<<showerSim.time()<<std::endl;
     }

     for (const auto& iter : showerMap)  _hCaShowerEdep->Fill(iter.first,iter.second);
     for (const auto& iter : showerMap2) _hCaShowerNum->Fill(iter.first,iter.second);


     //look at reconstructed hits
     if (!CaloHitsHandle.isValid()) return;

     double totalEdep = 0.0;
     set<int> hit_crystals;

     for (unsigned int ic=0; ic<CaloHits.size();++ic)
     {
         const CaloHit &hit     = CaloHits.at(ic);

         totalEdep += hit.energyDep();
         hit_crystals.insert(hit.crystalID());

	 if ( _diagLevel > 1 && _nAnalyzed < _maxFullPrint )
	   cout<<"Readback: caloHit id = "<<hit.crystalID()<<" "<<"energy = "<<hit.energyDep()<<" time= "<<hit.time()<<endl;
     }

     _hCaEdep->Fill(totalEdep);
     _hCaNcrystal->Fill(hit_crystals.size());


  }

  void ReadBack::doTracker(const art::Event& event){

    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = *GeomHandle<Tracker>();

    // Get a ValidHandle to the generated particles.
    auto gens = event.getValidHandle<GenParticleCollection>(_gensTag);

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel,simParticles);

    // Handle to information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoMultiCollection> volsHandle;
    event.getSubRun().getByLabel(_g4ModuleLabel,volsHandle);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volsHandle.isValid() );

    PhysicalVolumeInfoSingleStage const* vols = (volsHandle.isValid() && !volsHandle->empty()) ?
      &volsHandle->at(0) : nullptr;

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volsHandle->empty());
    }

    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hits->size());

    // A silly example just to show that we have a messsage logger.
    if ( hits->size() > 300 ){
      mf::LogWarning("HitInfo")
        << "Number of hits "
        << hits->size()
        << " may be too large.";
    }

    // A silly example just to show how to throw.
    if ( hits->size() > 1000000 ){
      throw cet::exception("RANGE")
        << "Way too many hits in this event.  Something is really wrong."
        << hits->size();
    }

    if ( _nAnalyzed < _maxFullPrint ){
      int ngen(0);
      for ( auto const& gen : *gens ){
        cout << "Readback Gen: "
             << event.id().event() << " "
             << ngen++ << " "
             << gen
             << endl;
      }
    }

    // ntuple buffer.
    float nt[_ntup->GetNvar()];

    // Loop over all hits.
    for ( size_t i=0; i<hits->size(); ++i ){

      // Alias, used for readability.
      const StepPointMC& hit = (*hits)[i];

      // Skip hits with low pulse height.
      if ( hit.eDep() < _minimumEnergy ) continue;

      // Get the hit information.
      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

      // Get the straw information:
      const Straw&      straw = tracker.getStraw( hit.strawId() );
      const CLHEP::Hep3Vector& mid   = straw.getMidPoint();
      const CLHEP::Hep3Vector& w     = straw.getDirection();

      // Count how many nearest neighbours are also hit.
      int nNeighbours = countHitNeighbours( straw, hits );

      // Compute an estimate of the drift distance.
      TwoLinePCA pca( mid, w, pos, mom);

      // Check that the radius of the reference point in the local
      // coordinates of the straw.  Should be 2.5 mm.
      double s = w.dot(pos-mid);
      CLHEP::Hep3Vector point = pos - (mid + s*w);

      // The simulated particle that made this hit.
      int trackId = hit.simParticle().key();

      double normPointMag = point.mag()/tracker.strawInnerRadius();
      double normS = s/straw.halfLength();

      if ( _diagLevel > 1 ){
        cout << __func__
             << " normalized reference point - 1 : "
             << scientific
             << normPointMag - 1.
             << " normalized wire z of reference point - 1 : "
             << std::abs(normS) -1.
             << fixed
             << endl;
      }

      if ( ( normPointMag - 1. > _strawHitPositionTolerance ) ||
           ( std::abs(normS) - 1. > _strawHitPositionTolerance ) ) {
        throw cet::exception("GEOM") << __func__
                                     << " Hit " << pos
                                     << " ouside the straw " << straw.id()
                                     << " inconsistent tracker geometry file? "
                                     << "; radial difference: "
                                     << point.mag()/tracker.strawInnerRadius() - 1.
                                     << ", longitudinal difference: "
                                     << std::abs(s)/straw.halfLength() - 1.
                                     << "; tolerance : "
                                     << _strawHitPositionTolerance
                                     << endl;
      }

      // Debug info

      //       // remove this for production, intended for transportOnly.py
      //       if (pca.dca()>strawDetail.innerRadius() || abs(point.mag()- strawDetail.innerRadius())>1.e-6 ) {

      //         cerr << "*** Bad hit?: "
      //              << event.id().event() << " "
      //              << i                  << " "
      //              << trackId            << "   "
      //              << hit.volumeId()     << " "
      //              << straw.id()         << " | "
      //              << pca.dca()          << " "
      //              << pos                << " "
      //              << mom                << " "
      //              << point.mag()        << " "
      //              << hit.eDep()         << " "
      //              << s
      //              << endl;
      //       }

      // Default values for these, in case information is not available.
      int pdgId(0);
      GenId genId;

      if ( haveSimPart ){
        SimParticle const& sim = *hit.simParticle();

        // PDG Particle Id of the sim particle that made this hit.
        pdgId = sim.pdgId();

        // If this is a generated particle, which generator did it come from?
        // This default constructs to "unknown".
        if ( sim.fromGenerator() ){
          GenParticle const& gen = *sim.genParticle();
          genId = gen.generatorId();
        }
      }

      // Fill some histograms
      _hRadius->Fill(pos.perp());
      _hEnergyDep->Fill(hit.eDep()/keV);
      _hTime->Fill(hit.time());
      _hHitNeighbours->Fill(nNeighbours);
      _hCheckPointRadius->Fill(point.mag());
      _hCheckPointRadiusW->Fill(normPointMag);
      _hCheckPointWireZ->Fill(normS);

      _hxHit->Fill(pos.x());
      _hyHit->Fill(pos.y());
      _hzHit->Fill(pos.z());

      _hDriftDist->Fill(pca.dca());
      _hDriftDistW->Fill(pca.dca()/tracker.strawInnerRadius());

      _hStepLength->Fill( hit.stepLength() );

      // Fill the ntuple.
      nt[0]  = event.id().event();
      nt[1]  = trackId;
      nt[2]  = hit.volumeId();
      nt[3]  = pos.x();
      nt[4]  = pos.y();
      nt[5]  = pos.z();
      nt[6]  = mid.x();
      nt[7]  = mid.y();
      nt[8]  = mid.z();
      nt[9]  = pca.dca();
      nt[10] = hit.time();
      nt[11] = straw.id().getPlane();
      nt[12] = straw.id().getPanel();
      nt[13] = straw.id().getLayer();
      nt[14] = pdgId;
      nt[15] = genId.id();
      nt[16] = hit.eDep()/keV;
      nt[17] = mom.mag();
      nt[18] = hit.stepLength();
      nt[19] = normS;
      nt[20] = straw.id().getStraw();

      _ntup->Fill(nt);

      // Fill the TGraph; need to manage size by hand.
      if ( _xyHitCount < _xyHitsMax ){
        _xyHits->SetPoint( _xyHitCount++, pos.x(), pos.y() );
      }

      // Print out limited to the first few events.
      if ( _nAnalyzed < _maxFullPrint ){

        cerr << "Readback"
             << " hit: "
             << event.id().event() << " "
             << i                  <<  " "
             << trackId               << "   "
             << straw.id().asUint16() << " "
             << straw.id()         << " | "
             << pca.dca()          << " "
             << pos                << " "
             << mom                << " "
             << point.mag()        << " "
             << hit.eDep()         << " "
             << hit.visibleEDep()  << " "
             << s
             << endl;
      }

    } // end loop over hits.

    // Additional printout and histograms about the simulated particles.
    if ( haveSimPart && (_nAnalyzed < _maxFullPrint) ){

      GlobalConstantsHandle<ParticleDataTable> pdt;

      for ( SimParticleCollection::const_iterator i=simParticles->begin();
            i!=simParticles->end(); ++i ){

        SimParticle const& sim = i->second;

        if ( sim.madeInG4() ) {

          _hMomentumG4->Fill( sim.startMomentum().rho() );

        } else {

          // Particle Data Group Id number of this SimParticle
          int pdgId = sim.pdgId();

          // Name of this particle type.
          ParticleDataTable::maybe_ref particle = pdt->particle(pdgId);
          string pname = particle.ref().name();

          // Information about generated particle.
          GenParticle const& gen = *sim.genParticle();
          GenId genId(gen.generatorId());

          // Physical volume in which this track started.
          PhysicalVolumeInfo const& volInfo = vols->at(cet::map_vector_key(sim.startVolumeIndex()));
          PhysicalVolumeInfo const& endInfo = vols->at(cet::map_vector_key(sim.endVolumeIndex()));

          cerr << "Readback"
               << " Simulated Particle: "
               << i->first            << " "
               << pdgId               << " "
               << genId.name()        << " "
               << sim.startPosition() << " "
               << volInfo.name()      << " "
               << volInfo.copyNo()    << " | "
               << sim.endPosition()   << " "
               << sim.endGlobalTime() << " "
               << endInfo.name()      << " "
               << endInfo.copyNo()    << " | "
               << sim.stoppingCode()  << " "
               << sim.startMomentum().vect().mag()
               << endl;
        }

      }
    }

  } // end doTracker

  // Count how many of this straw's nearest neighbours are hit.
  // If we have enough hits per event, it will make sense to make
  // a class to let us direct index into a list of which straws have hits.
  int ReadBack::countHitNeighbours( Straw const& straw,
                                    art::Handle<StepPointMCCollection>& hits ){

    int count(0);

    for ( auto const& step : *hits ){
      if ( step.strawId().nearestNeighbor(straw.id()) ){
	++count;
      }
    }
    return count;
  }  // end countHitNeighbours

  //
  // Example of how to read information about stopping target
  //
  // Here we assume that the primary particles are the conversion
  // electrons generated in the stopping target.
  //
  void ReadBack::doStoppingTarget(const art::Event& event) {

    // Find original G4 steps in the stopping target
    art::Handle<StepPointMCCollection> sthits;
    event.getByLabel(_g4ModuleLabel,_targetStepPoints,sthits);

    // SimParticles container
    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel,simParticles);
    if( !(simParticles.isValid()) || simParticles->empty() ) return;

    // Loop over all hits in the stopping target. Check, that the
    // hit belongs to primary particle. If so, calculate total energy
    // deposition in the target, total path length and the number of
    // foils, crossed by the electron.

    int id_start = -1;
    double eDep=0.0, pathLength=0.0;
    map<int,int> foils;

    // Loop over all hits in the stopping target
    for ( size_t i=0; i<sthits->size(); ++i ){

      // This is G4 hit (step) in the target
      const StepPointMC& hit = (*sthits)[i];

      // Here we select only those hits, which are generated by
      // the primary track - which is assumed to be the original
      // particle, generated by ConversionGun
      SimParticleCollection::key_type trackId = hit.trackId();
      if( trackId.asInt() != 1 ) continue;

      // Here we require that there is information about the primary
      // particle in the SimParticle collection. It is not neccessary for
      // this example, but it is typical requirement in the real analysis
      SimParticle const* sim = simParticles->getOrNull(trackId);
      if( !sim ) continue;

      // Get the foil id where the hit occured. If it is the first hit,
      // remember this id as the source foil id.
      int id = hit.volumeId();
      if( id_start<0 ) id_start=id;

      // Here we calculate number of steps in each foil. There could be
      // many hits in each foil. This number is not used in the example,
      // it is calculated here just as an example. But we use this map to
      // calculate number of foils with the hits. Be aware, that if particle
      // crosses foil without energy deposition (just passes it), it still
      // counts as a hit. Here we record all foils particle crosses.
      // If we want to record only those foils where particle had
      // interactions, we would need to do the following:
      //      if( hit.totalEDep()>0 ) foils[id]++;
      foils[id]++;

      // Calculate total energy loss and path length primary particle
      // has in the target.
      eDep += hit.totalEDep();
      pathLength += hit.stepLength();

    }

    // Number of crossed foils
    int nfoil = foils.size();

    // Fill histograms
    if( id_start>=0 ) {
      _hTargetEdep->Fill(eDep);
      _hTargetPathLength->Fill(pathLength);
      _hTargetNfoils->Fill(nfoil);
      _hTargetNfoils2D->Fill(id_start,nfoil);
    }

  } // end doStoppingTarget

  void ReadBack::endJob(){
    cout << "ReadBack::endJob Number of events skipped "
         << "due to G4 completion status: "
         << _nBadG4Status
         << endl;
  }

  void ReadBack::doCRV(const art::Event& event){

    // Get a reference to CosmicRayShield (it contains crv)

    GeomHandle<CosmicRayShield> cosmicRayShieldGeomHandle;
    CosmicRayShield const& crv(*cosmicRayShieldGeomHandle);

    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( _g4ModuleLabel, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_crvStepPoints,hits);
    if ( ! hits.isValid() ) { return; }

    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> genParticles;
    event.getByLabel(_generatorModuleLabel,genParticles);

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel,simParticles);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() );
    }

    // A silly example just to show that we have a messsage logger.
    if ( hits->size() > 300 ){
      mf::LogWarning("HitInfo")
        << "Number of CRV hits "
        << hits->size()
        << " may be too large.";
    }

    // A silly example just to show how to throw.
    if ( hits->size() > 1000000 ){
      throw cet::exception("RANGE")
        << "Way too many CRV hits in this event.  Something is really wrong."
        << hits->size();
    }

    if ( _nAnalyzed < _maxFullPrint ){
      cerr << g4Status << endl;
    }

    // ntuple buffer.
    float nt[_ntupCRV->GetNvar()];

    // Loop over all hits.
    for ( size_t i=0; i<hits->size(); ++i ){

      // Alias, used for readability.
      const StepPointMC& hit = (*hits)[i];

      // Get the hit information.
      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

      // Get the CRSScintillatorBar information:
      const CRSScintillatorBar&  bar = crv.getBar( hit.barIndex() );
      CLHEP::Hep3Vector const &  mid = bar.getPosition();

      CLHEP::Hep3Vector barLengths = CLHEP::Hep3Vector(bar.getHalfLengths()[0],
                                                       bar.getHalfLengths()[1],
                                                       bar.getHalfLengths()[2]);

      CLHEP::Hep3Vector hitLocal  = pos-mid;
      CLHEP::Hep3Vector hitLocalN = CLHEP::Hep3Vector(hitLocal.x()/barLengths.x(),
                                                      hitLocal.y()/barLengths.y(),
                                                      hitLocal.z()/barLengths.z());

      // The simulated particle that made this hit.
      SimParticleCollection::key_type trackId(hit.trackId());

      // Default values for these, in case information is not available.
      int pdgId(0);
      GenId genId;

      if ( haveSimPart ){

        SimParticle const& sim = simParticles->at(trackId);

        // PDG Particle Id of the sim particle that made this hit.
        pdgId = sim.pdgId();

        // If this is a generated particle, which generator did it come from?
        // This default constructs to "unknown".
        if ( sim.fromGenerator() ){
          GenParticle const& gen = *sim.genParticle();
          genId = gen.generatorId();
        }

      }

      // Fill the ntuple.
      nt[ 0] = event.id().event();
      nt[ 1] = hit.trackId().asInt();
      nt[ 2] = hit.volumeId();
      nt[ 3] = pos.x();
      nt[ 4] = pos.y();
      nt[ 5] = pos.z();
      nt[ 6] = mid.x();
      nt[ 7] = mid.y();
      nt[ 8] = mid.z();
      nt[ 9] = fabs((pos-mid).x());
      nt[10] = fabs((pos-mid).y());
      nt[11] = fabs((pos-mid).z());
      nt[12] = hitLocalN.x();
      nt[13] = hitLocalN.y();
      nt[14] = hitLocalN.z();
      nt[15] = hit.time();
      nt[16] = bar.id().getShieldNumber();
      nt[17] = bar.id().getModuleNumber();
      nt[18] = bar.id().getLayerNumber();
      nt[19] = pdgId;
      nt[20] = genId.id();
      nt[21] = hit.eDep()/keV;
      nt[22] = mom.mag();
      nt[23] = hit.stepLength();

      _ntupCRV->Fill(nt);

      _hCRVMultiplicity->Fill(hit.volumeId());

      // Print out limited to the first few events.
      if ( _nAnalyzed < _maxFullPrint ){

        cerr << "Readback"
             << " hit: "
             << event.id().event() << " "
             << i                  <<  " "
             << hit.trackId()      << "   "
             << hit.volumeId()     << " "
             << bar.id()           << " | "
             << pos                << " "
             << mid                << " "
             << (mid-pos)          << " | "
             << hitLocalN          << " | "
             << mom                << " "
             << hit.totalEDep()/keV << " "
             << hit.stepLength()
             << endl;
      }

    } // end loop over hits.

  } // end doCRV

  void ReadBack::doPointTrajectories ( const art::Event& event ){

    art::Handle<MCTrajectoryCollection> trajsHandle;
    event.getByLabel(_g4ModuleLabel,trajsHandle);
    MCTrajectoryCollection const& trajs(*trajsHandle);

    _hNPointTrajectories->Fill( trajs.size() );

    for ( MCTrajectoryCollection::const_iterator i=trajs.begin();
          i != trajs.end(); ++i ){
      const MCTrajectory& traj(i->second);

      _hNPointsPerTrajectory->Fill( traj.size() );

      double logSize = ( traj.size() == 0 ) ? -1 : log10(traj.size());
      _hNPointsPerTrajectoryLog->Fill( logSize );
    }

  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ReadBack);
