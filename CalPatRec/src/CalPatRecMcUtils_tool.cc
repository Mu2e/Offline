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

#include "CalPatRec/inc/McUtilsToolBase.hh"

namespace mu2e {

  class CalPatRecMcUtils : public mu2e::McUtilsToolBase {
  public:

    CalPatRecMcUtils(const fhicl::ParameterSet& PSet);
    ~CalPatRecMcUtils();

  public:

    virtual double mcDoca(const art::Event* Event     , 
			  const char*       MCCollName, 
			  const Straw*      Straw     ) override ;

    virtual int    nGenHits(const art::Event*         Event         , 
			    fhicl::ParameterSet*      TimeOffsets   ,
			    const char*               MCDigiCollName, 
			    const StrawHitCollection* Shcol         ) override;

    virtual const StrawDigiMCCollection* getListOfMcStrawHits(const art::Event* Event,
							      const art::InputTag& Tag) override;
    
    virtual const SimParticle* getSimParticle(const StrawDigiMCCollection* List, int IHit) override;

    int   getID      (const SimParticle* Sim) override;
    int   getPdgID   (const SimParticle* Sim) override;
    float getStartMom(const SimParticle* Sim) override;
  };

//-----------------------------------------------------------------------------
  CalPatRecMcUtils::CalPatRecMcUtils(const fhicl::ParameterSet& PSet) {
  }

//-----------------------------------------------------------------------------
  CalPatRecMcUtils::~CalPatRecMcUtils() {
  }


//-----------------------------------------------------------------------------
// find MC truth DOCA in a given straw
// start from finding the right vector of StepPointMC's
//-----------------------------------------------------------------------------
  double CalPatRecMcUtils::mcDoca(const art::Event* Event, const char* MCCollName, const Straw* Straw) {

    static int    last_event(-1);
    //    static int    first_call( 1);

    static const StrawDigiMCCollection*  listOfMCStrawHits(NULL);
    

    double mcdoca(-99.0);

    int iev = Event->event();

    if (iev != last_event) {
      art::Handle<mu2e::StrawDigiMCCollection> mcdigiH;
      Event->getByLabel(MCCollName,mcdigiH);
      if (mcdigiH.isValid()) listOfMCStrawHits = (mu2e::StrawDigiMCCollection*) mcdigiH.product();
      else                   listOfMCStrawHits = NULL;

      last_event = iev;
    }

    if (listOfMCStrawHits) { 
      int nstraws = listOfMCStrawHits->size();

      const mu2e::StepPointMC* step(0);

      for (int i=0; i<nstraws; i++) {
	const mu2e::StrawDigiMC*  mcdigi = &listOfMCStrawHits->at(i);

	if (mcdigi->wireEndTime(mu2e::TrkTypes::cal) < mcdigi->wireEndTime(mu2e::TrkTypes::hv)) {
	  step = mcdigi->stepPointMC(mu2e::TrkTypes::cal).get();
	}
	else {
	  step = mcdigi->stepPointMC(mu2e::TrkTypes::hv ).get();
	}

	int volume_id = step->volumeId();
	if (volume_id == Straw->id().asUint16()) {
//-----------------------------------------------------------------------------
// step found - use the first one in the straw
//-----------------------------------------------------------------------------
	  break;
	}
      }

      if (step) {
	const CLHEP::Hep3Vector* v1 = &Straw->getMidPoint();
	HepPoint p1(v1->x(),v1->y(),v1->z());

	const CLHEP::Hep3Vector* v2 = &step->position();
	HepPoint    p2(v2->x(),v2->y(),v2->z());

	TrkLineTraj trstraw(p1,Straw->getDirection()  ,0.,0.);
	TrkLineTraj trstep (p2,step->momentum().unit(),0.,0.);

	TrkPoca poca(trstep, 0., trstraw, 0.);

	mcdoca = poca.doca();
      }
    }

    return mcdoca;
  }

//-----------------------------------------------------------------------------
// calculates N(MC hits) produced by the signal particle, SIM_ID = 1, with P > 100
//-----------------------------------------------------------------------------
  int CalPatRecMcUtils::nGenHits(const art::Event*         Event         ,
				 fhicl::ParameterSet*      TimeOffsets   ,
				 const char*               MCDigiCollName,
				 const StrawHitCollection* Shcol         ) {

    static int     last_event(-1);
    static int     first_call(1);
    static double  mbtime;

    static SimParticleTimeOffset*                 timeOffsets(NULL);
    static const PtrStepPointMCVectorCollection*  listOfMCStrawHits(NULL);

    double  time_threshold(500.);
    int     n_gen_hits(  0 );

    if (first_call == 1) {
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime = accPar->deBuncherPeriod;

      timeOffsets = new SimParticleTimeOffset(*TimeOffsets);
      first_call  = 0;
    }
//-----------------------------------------------------------------------------
// update if new event
//-----------------------------------------------------------------------------
    int iev = Event->event();

    if (iev != last_event) {
      art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
      Event->getByLabel(MCDigiCollName,mcptrHandle);
      if (mcptrHandle.isValid()) listOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
      else                       listOfMCStrawHits = NULL;

      timeOffsets->updateMap(*Event);

      last_event = iev;
    }

    if (listOfMCStrawHits == NULL) return -1;

    double  pEntrance(.0);

    int nhits = Shcol->size();
    for (int i=0; i<nhits; i++) {

      const mu2e::StrawDigiMC* mcdigi = &listOfMCStrawHits->at(i);

      const mu2e::StepPointMC   *step;
      if (mcdigi->wireEndTime(mu2e::TrkTypes::cal) < mcdigi->wireEndTime(mu2e::TrkTypes::hv)) {
	step = mcdigi->stepPointMC(mu2e::TrkTypes::cal).get();
      }
      else {
	step = mcdigi->stepPointMC(mu2e::TrkTypes::hv ).get();
      }

      int gen_index(-1), sim_id(-1);

      if (step) {
	art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle();

	if (simptr->fromGenerator()) gen_index = simptr->genParticle()->generatorId().id();
	else                         gen_index = -1;

	sim_id        = simptr->id().asInt();
      }

      if ((gen_index > 0) && (sim_id == 1)) {
	double step_time = timeOffsets->timeWithOffsetsApplied(*step);
	step_time = fmod(step_time,mbtime);
	if (step_time > time_threshold) {
	  ++n_gen_hits;
	  double pstep = step->momentum().mag();
	  if (pstep > pEntrance) {
	    pEntrance = pstep;
	  }
	}
      }
    }

    if (pEntrance < 100. ) n_gen_hits = 0;

    return n_gen_hits;
  }

//-----------------------------------------------------------------------------
  const StrawDigiMCCollection* CalPatRecMcUtils::getListOfMcStrawHits(const art::Event* Event,const art::InputTag& Tag) {
    auto handle = Event->getValidHandle<StrawDigiMCCollection>(Tag);
    const StrawDigiMCCollection* coll = handle.product();
    return coll;
  }

//-----------------------------------------------------------------------------
  const SimParticle* CalPatRecMcUtils::getSimParticle(const StrawDigiMCCollection* List, int IHit) {
    const mu2e::StrawDigiMC* mcdigi = &List->at(IHit);

    const mu2e::StepPointMC   *step;
    if (mcdigi->wireEndTime(mu2e::TrkTypes::cal) < mcdigi->wireEndTime(mu2e::TrkTypes::hv)) {
      step = mcdigi->stepPointMC(mu2e::TrkTypes::cal).get();
    }
    else {
      step = mcdigi->stepPointMC(mu2e::TrkTypes::hv ).get();
    }

    const mu2e::SimParticle* sim = &(*step->simParticle());

    return sim;
  }

//-----------------------------------------------------------------------------
  int   CalPatRecMcUtils::getID      (const SimParticle* Sim) { return Sim->id().asInt();  }

//-----------------------------------------------------------------------------
  int   CalPatRecMcUtils::getPdgID   (const SimParticle* Sim) { return Sim->pdgId();  }

//-----------------------------------------------------------------------------
  float CalPatRecMcUtils::getStartMom(const SimParticle* Sim) {
    CLHEP::HepLorentzVector const& p = Sim->startMomentum();
    return sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());
  }

  DEFINE_ART_CLASS_TOOL(CalPatRecMcUtils)
}
