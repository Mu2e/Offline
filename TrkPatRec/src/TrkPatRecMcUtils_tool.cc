///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"

#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/BbrGeom/TrkLineTraj.hh"

#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "CLHEP/Matrix/Vector.h"

#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"

#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"

#include "TrackerGeom/inc/Straw.hh"

#include "Mu2eUtilities/inc/McUtilsToolBase.hh"

namespace mu2e {

  class TrkPatRecMcUtils : public mu2e::McUtilsToolBase {
  protected:
    std::string                   _strawDigiMCCollTag;
    const StrawDigiMCCollection*  _mcdigis;
    unsigned int                  _lastEvent;
    //    SimParticleTimeOffset*       _timeOffsets;
    double                        _mbtime;

  public:

    TrkPatRecMcUtils(const fhicl::ParameterSet& PSet);
    ~TrkPatRecMcUtils();

  public:

    int     initEvent(const art::Event* Event);

    virtual int    strawHitSimId(const art::Event* Event, int Index) override;

    virtual double mcDoca(const art::Event* Event, 
			  int               Index, 
			  const Straw*      Straw) override ;

    // virtual int    nGenHits(const art::Event*         Event      , 
    // 			    fhicl::ParameterSet*      TimeOffsets,
    // 			    const StrawHitCollection* Shcol      ) override;

    virtual const StrawDigiMCCollection* getListOfMcStrawHits(const art::Event* Event,
							      const art::InputTag& Tag) override;
    
    virtual const SimParticle* getSimParticle(const art::Event* Event, int HitIndex) override;

    int   getID      (const SimParticle* Sim) override;
    int   getPdgID   (const SimParticle* Sim) override;
    float getStartMom(const SimParticle* Sim) override;

    const StepPointMC* getStepPointMC(int HitIndex);
  };

//-----------------------------------------------------------------------------
  TrkPatRecMcUtils::TrkPatRecMcUtils(const fhicl::ParameterSet& PSet) : 
    _strawDigiMCCollTag{ PSet.get<std::string>("strawDigiMCCollTag") }
  {
    _lastEvent   = -1;
    _mcdigis     = nullptr;
    _mbtime      = -1;
    //    _timeOffsets = new SimParticleTimeOffset(*TimeOffsets);
  }

//-----------------------------------------------------------------------------
  TrkPatRecMcUtils::~TrkPatRecMcUtils() {
    //    delete _timeOffsets;
  }

//-----------------------------------------------------------------------------
  int TrkPatRecMcUtils::initEvent(const art::Event* Event) {
    art::Handle<mu2e::StrawDigiMCCollection> mcdigiH;
    Event->getByLabel(_strawDigiMCCollTag,mcdigiH);
    if (mcdigiH.isValid()) _mcdigis = (mu2e::StrawDigiMCCollection*) mcdigiH.product();
    else                   _mcdigis = NULL;

    //    _timeOffsets->updateMap(*Event);
    if (_mbtime < 0) {
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      _mbtime      = accPar->deBuncherPeriod;
    }
    _lastEvent = Event->event();

    return 0;
  }

//-----------------------------------------------------------------------------
  const StepPointMC* TrkPatRecMcUtils::getStepPointMC(int HitIndex) {
    const StepPointMC* step(nullptr);
    if (_mcdigis) { 
      const mu2e::StrawDigiMC*  mcdigi = &_mcdigis->at(HitIndex);
      
      if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
	step = mcdigi->stepPointMC(mu2e::StrawEnd::cal).get();
      }
      else {
	step = mcdigi->stepPointMC(mu2e::StrawEnd::hv ).get();
      }
    }
    return step;
  }

//-----------------------------------------------------------------------------
// returns ID of the SimParticle corresponding to straw hit 'Index'
//-----------------------------------------------------------------------------
  int TrkPatRecMcUtils::strawHitSimId(const art::Event* Event, int HitIndex) {
    int           simId(-1);
    
    if (Event->event() != _lastEvent) initEvent(Event);
    
    const StepPointMC* step = getStepPointMC(HitIndex);

    if (step) simId = step->simParticle()->id().asInt();
    
    return simId;
  }
//-----------------------------------------------------------------------------
// find MC truth DOCA in a given straw
// start from finding the right vector of StepPointMC's
//-----------------------------------------------------------------------------
  double TrkPatRecMcUtils::mcDoca(const art::Event* Event, int HitIndex, const Straw* Straw) {

    double mcdoca(-99.0);

    if (Event->event() != _lastEvent) initEvent(Event);

    const CLHEP::Hep3Vector* v1 = &Straw->getMidPoint();
    HepPoint p1(v1->x(),v1->y(),v1->z());

    const StepPointMC* step = getStepPointMC(HitIndex);

    const CLHEP::Hep3Vector* v2 = &step->position();
    HepPoint    p2(v2->x(),v2->y(),v2->z());

    TrkLineTraj trstraw(p1,Straw->getDirection()  ,0.,0.);
    TrkLineTraj trstep (p2,step->momentum().unit(),0.,0.);

    TrkPoca poca(trstep, 0., trstraw, 0.);

    mcdoca = poca.doca();

    return mcdoca;
  }

// //-----------------------------------------------------------------------------
// // calculates N(MC hits) produced by the signal particle, SIM_ID = 1, with P > 100
// //-----------------------------------------------------------------------------
//   int TrkPatRecMcUtils::nGenHits(const art::Event*         Event         ,
// 				 fhicl::ParameterSet*      TimeOffsets   ,
// 				 const StrawHitCollection* Shcol         ) {

//     //    static int     last_event(-1);
//     //    static int     first_call(1);

//     double  time_threshold(500.);
//     int     n_gen_hits(  0 );
// //-----------------------------------------------------------------------------
// // update if new event
// //-----------------------------------------------------------------------------
//     if (Event->event() != _lastEvent) initEvent(Event);

//     if (_mcdigis == NULL) return -1;

//     double  pEntrance(.0);

//     int nhits = Shcol->size();
//     for (int i=0; i<nhits; i++) {
//       const mu2e::StepPointMC*   *step = getStepPointMC(i);

//       int gen_index(-1), sim_id(-1);

//       if (step) {
// 	art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle();

// 	if (simptr->fromGenerator()) gen_index = simptr->genParticle()->generatorId().id();
// 	else                         gen_index = -1;

// 	sim_id        = simptr->id().asInt();
//       }

//       if ((gen_index > 0) && (sim_id == 1)) {
// 	double step_time = timeOffsets->timeWithOffsetsApplied(*step);
// 	step_time = fmod(step_time,_mbtime);
// 	if (step_time > time_threshold) {
// 	  ++n_gen_hits;
// 	  double pstep = step->momentum().mag();
// 	  if (pstep > pEntrance) {
// 	    pEntrance = pstep;
// 	  }
// 	}
//       }
//     }

//     if (pEntrance < 100. ) n_gen_hits = 0;

//     return n_gen_hits;
//   }

//-----------------------------------------------------------------------------
  const StrawDigiMCCollection* TrkPatRecMcUtils::getListOfMcStrawHits(const art::Event* Event,const art::InputTag& Tag) {
    auto handle = Event->getValidHandle<StrawDigiMCCollection>(Tag);
    const StrawDigiMCCollection* coll = handle.product();
    return coll;
  }

//-----------------------------------------------------------------------------
  const SimParticle* TrkPatRecMcUtils::getSimParticle(const art::Event* Event, int HitIndex) {

    if (Event->event() != _lastEvent) initEvent(Event);

    const mu2e::StepPointMC   *step = getStepPointMC(HitIndex);

    const mu2e::SimParticle* sim(nullptr);

    sim = (step != nullptr) ? step->simParticle().get() : nullptr;

    return sim;
  }

//-----------------------------------------------------------------------------
  int   TrkPatRecMcUtils::getID      (const SimParticle* Sim) { return Sim->id().asInt();  }

//-----------------------------------------------------------------------------
  int   TrkPatRecMcUtils::getPdgID   (const SimParticle* Sim) { return Sim->pdgId();  }

//-----------------------------------------------------------------------------
  float TrkPatRecMcUtils::getStartMom(const SimParticle* Sim) {
    CLHEP::HepLorentzVector const& p = Sim->startMomentum();
    return sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());
  }

  DEFINE_ART_CLASS_TOOL(TrkPatRecMcUtils)
}
