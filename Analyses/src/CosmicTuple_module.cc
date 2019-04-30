//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// $Id: CosmicTuple_module.cc,v 1.11 2013/09/27 16:03:41 gandr Exp $
// $Author: gandr $
// $Date: 2013/09/27 16:03:41 $
//
// Original author Yury Kolomensky (Rob Kutschke)
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>
#include <string>
#include <vector>

using namespace std;

namespace mu2e {

  class CosmicTuple : public art::EDFilter {
  public:

    explicit CosmicTuple(fhicl::ParameterSet const& pset);
    virtual ~CosmicTuple() { }

    virtual void beginJob();
    virtual bool beginRun(art::Run &r);

    // This is called for each event.
    virtual bool filter(art::Event& e);

  private:

    // Module label of the geerator module.
    std::string _generatorModuleLabel;

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

  CosmicTuple::CosmicTuple(fhicl::ParameterSet const& pset) :
    EDFilter{pset},
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _minimump(pset.get<double>("minimump")),
    _maximump(pset.get<double>("maximump")),
    _minHits(pset.get<int>("minHits")),

    _nAnalyzed(0),
    _hEventsize(0),
    _ntupTrk(0)

  {
  }

  void CosmicTuple::beginJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    // Create some 1D histograms.
    _hEventsize     = tfs->make<TH1F>( "EventSize",       "Size of the Event",     100,  0., 100000. );

    // Create an ntuple.
    _ntupTrk = tfs->make<TNtuple>( "ntupTrk", "Trk ntuple",
"evt:trk:pid:pmag:genId:pidGen:eGen:thGen:xGen:yGen:zGen:nHits:ptrs:xtrs:ytrs:ztrs:pprs:xprs:yprs:zprs:prnId:time:calE:nCryst:nAPD:pAng:prCrea:prStop:trCrea:trStop:run:px:py:pz:E:vx:vy:vz:vt:isSh:ppre:xpre:ypre:zpre:trvs:trve:prvs:prve:calEi");
  }

  bool CosmicTuple:: beginRun(art::Run& run) {
    _runNumber = run.id().run();
    return true;
  }

  bool CosmicTuple::filter(art::Event& event) {

    typedef SimParticleCollection::key_type key_type;

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    if ( _nAnalyzed % 10000 == 0 ) {
      mf::LogInfo("CosmicTuple")
        << "Processing event " << _nAnalyzed;
    }

    // Ask the event to give us a "handle" to the requested hits.
       static const string collectionName("tracker");
       art::Handle<StepPointMCCollection> hits;
       event.getByLabel(_g4ModuleLabel,collectionName,hits);

    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> genParticles;
    event.getByLabel(_generatorModuleLabel, genParticles);

    // Get handles to the generated and simulated particles.
    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() );
    }

    // ntuple buffer.
    //float ntT[_ntupTrk->GetNvar()];

    bool    pass     = false;
    /*
    bool    isSh     = false;
    double  calEne   = 0.;
    double  calEind  = 0.;
    int     prntPdg  = 0;
    int     prCr     = 0;
    int     prSt     = 0;
    int     trCr     = 0;
    int     trSt     = 0;
    double  px       = 0.;
    double  py       = 0.;
    double  pz       = 0.;
    double  E        = 0.;
    double  vx       = 0.;
    double  vy       = 0.;
    double  vz       = 0.;
    double  vt       = 0.;
    double  rmass    = 0.;
    unsigned trSVolume = 0;
    unsigned trEVolume = 0;
    unsigned prSVolume = 0;
    unsigned prEVolume = 0;
    */

    map<int,int> hit_crystals;
    map<int,int> hit_apds;

    art::ServiceHandle<GeometryService> geom;
    return pass;
    /*
    GeomHandle<VaneCalorimeter> cg;

      // Fill some histograms
      _hEventsize->Fill(simParticles->size());

    if ( simParticles->size() >= 50000 ) return pass;

    art::Handle<StepPointMCCollection> rohits;
    event.getByLabel(_g4ModuleLabel,"calorimeter",rohits);

    art::Handle<StepPointMCCollection> apdhits;
    event.getByLabel(_g4ModuleLabel,"calorimeterRO",apdhits);

    GlobalConstantsHandle<ParticleDataTable> pdt;

    // Construct an object that ties together all of the simulated particle and hit info.
    SimParticlesWithHits sims( event,
                               "g4run",
                               "makeSH",
                               "tracker",
                               0.001,
                               0 );

    typedef SimParticlesWithHits::map_type map_type;

    for ( map_type::const_iterator i=sims.begin();
          i != sims.end(); ++i ){

      // All information about this SimParticle
      SimParticleInfo const& simInfo = i->second;
      SimParticle const& sim = simInfo.simParticle();

      // Information about primary ancestor
      SimParticleAncestors ancestor( sim,
                                    *simParticles,
                                    *genParticles);
      GenParticle const& gen_parent = ancestor.originalGen();

      // Immediate parent particle
      SimParticle const* sim_parent = 0;
      if( sim.hasParent() ) sim_parent = simParticles->getOrNull(sim.parentId());

      // Information about StrawHits that belong on this SimParticle.
      vector<StrawHitMCInfo> const& infos = simInfo.strawHitInfos();
      int nHits = infos.size();

      CLHEP::Hep3Vector y(0,-1,0);
      CLHEP::Hep3Vector posGen = gen_parent.position();
      CLHEP::Hep3Vector trspos(0,0,0);
      CLHEP::Hep3Vector prspos(0,0,0);
      CLHEP::Hep3Vector trsmom(0,0,0);
      CLHEP::Hep3Vector prsmom(0,0,0);
      CLHEP::Hep3Vector prepos(0,0,0);
      CLHEP::Hep3Vector premom(0,0,0);

      trspos = sim.startPosition();
      trsmom = sim.startMomentum();
      trSVolume = sim.startVolumeIndex();
      trEVolume = sim.endVolumeIndex();
      ProcessCode creationCode = sim.creationCode();
      ProcessCode stoppingCode = sim.stoppingCode();
      trCr = creationCode;
      trSt = stoppingCode;

      if( sim.hasParent()) {
        prspos = sim_parent->startPosition();
        prsmom = sim_parent->startMomentum();
        prepos = sim_parent->endPosition();
        premom = sim_parent->endMomentum();
        prSVolume = sim_parent->startVolumeIndex();
        prEVolume = sim_parent->endVolumeIndex();
        prntPdg= sim_parent->pdgId();
        prCr   = sim_parent->creationCode();
        prSt   = sim_parent->stoppingCode();
      } else {
        prspos  = posGen;
        prsmom  = gen_parent.momentum();
        prntPdg = gen_parent.pdgId();
      }

      double  momentum = -1.;
      double  pitchAng = -1.;

     //Get particle mass
     ParticleDataTable::maybe_ref e_data = pdt->particle(sim.pdgId());
     if ( e_data ){
       rmass = e_data.ref().mass().value();
     } else{

       // If particle is unknown, set rest mass to 0.
       rmass = 0.;
       mf::LogWarning("PDGID")
         << "Particle ID created by G4 not known to the particle data table:  "
         << sim.pdgId()
         << "\n";
     }


      // First  StrawsHits to which this SimParticle contributed.
      StrawHitMCInfo const& info = infos.at(0);

      // Loop over all StepPointMC's that contribute to this StrawHit.
      std::vector<StepPointMC const *> const& steps = info.steps();

      for ( size_t k=0; k<steps.size(); ++k){
        StepPointMC const& step = *(steps[k]);

        if(step.trackId() == sim.id()) {
          momentum = step.momentum().mag();
          pitchAng = acos(step.momentum().z()/step.momentum().mag())*180./3.14159;
          px = step.momentum().x();
          py = step.momentum().y();
          pz = step.momentum().z();
          E  = sqrt(step.momentum().mag()*step.momentum().mag() + rmass*rmass);
          vx = step.position().x();
          vy = step.position().y();
          vz = step.position().z();
          vt = step.time();
          isSh = info.isShared();
          break;
        }
      }

      if( rohits.isValid() ) {
        calEne = 0.0;
        calEind= 0.0;
        for ( size_t i=0; i<rohits->size(); ++i ) {
          const StepPointMC & rohit = rohits->at(i);
          calEne += rohit.eDep();
          int cid = rohit.volumeId();
          hit_crystals[cid] =1;
          SimParticle const * csim = simParticles->getOrNull(rohit.trackId());
          while ( csim && csim->id() != sim.id() ) {
            csim = simParticles->getOrNull(csim->parentId());
          }
          if(csim){
            calEind += rohit.eDep();
          }
        }
      } else{
        calEne = -1;
        calEind= -1;
      }

      // Find original G4 steps in the APDs
      if( apdhits.isValid() ) {
        for ( size_t i=0; i<apdhits->size(); ++i ) {
          const StepPointMC & apdhit = apdhits->at(i);
          int apdid = apdhit.volumeId();
          int cida  = cg->caloInfo().crystalByRO(apdid);
          hit_apds[cida] =1;
        }
      }

      if( momentum < _minimump || momentum > _maximump ) continue;
      if( nHits< _minHits ) continue;
      pass = true;

      ntT[0]  = event.id().event();
      ntT[1]  = sim.id().asInt();
      ntT[2]  = sim.pdgId();
      ntT[3]  = momentum;
      ntT[4]  = gen_parent.generatorId().id();
      ntT[5]  = gen_parent.pdgId();
      ntT[6]  = gen_parent.momentum().e();
      ntT[7]  = y.angle(gen_parent.momentum().vect());
      ntT[8]  = posGen.x();
      ntT[9]  = posGen.y();
      ntT[10] = posGen.z();
      ntT[11] = nHits;
      ntT[12] = trsmom.mag();
      ntT[13] = trspos.x();
      ntT[14] = trspos.y();
      ntT[15] = trspos.z();
      ntT[16] = prsmom.mag();
      ntT[17] = prspos.x();
      ntT[18] = prspos.y();
      ntT[19] = prspos.z();
      ntT[20] = prntPdg;
      ntT[21] = sim.startGlobalTime();
      ntT[22] = calEne;
      ntT[23] = hit_crystals.size();
      ntT[24] = hit_apds.size();
      ntT[25] = pitchAng;
      ntT[26] = prCr;
      ntT[27] = prSt;
      ntT[28] = trCr;
      ntT[29] = trSt;
      ntT[30] = _runNumber;
      ntT[31] = px;
      ntT[32] = py;
      ntT[33] = pz;
      ntT[34] = E;
      ntT[35] = vx;
      ntT[36] = vy;
      ntT[37] = vz;
      ntT[38] = vt;
      ntT[39] = isSh;
      ntT[40] = premom.mag();
      ntT[41] = prepos.x();
      ntT[42] = prepos.y();
      ntT[43] = prepos.z();
      ntT[44] = trSVolume ;
      ntT[45] = trEVolume ;
      ntT[46] = prSVolume ;
      ntT[47] = prEVolume ;
      ntT[48] = calEind;

      _ntupTrk->Fill(ntT);
      

    } // end loop over hits.

    */
    return pass;

  } // end filter

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CosmicTuple);
