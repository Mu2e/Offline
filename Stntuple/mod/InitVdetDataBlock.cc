///////////////////////////////////////////////////////////////////////////////
// 2014-01-26 P.Murat
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/obj/TVdetDataBlock.hh"

//#include "RecoDataProducts/inc/StrawHitCollection.hh"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "Stntuple/mod/StntupleGlobals.hh"

mu2e::SimParticleTimeOffset* fgTimeOffsets;

//-----------------------------------------------------------------------------
Int_t StntupleInitMu2eVirtualDataBlock(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) 
{
  // initialize COT data block with the `event' data
  // return -1, if bank doesn't exist, 0, if everything is OK

  int ev_number, rn_number, nhits;

  static char   strh_module_label[100], strh_description[100];

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

  TVdetDataBlock* data = (TVdetDataBlock*) Block;
  data->Clear();
//-----------------------------------------------------------------------------
//  virtual hit information
//-----------------------------------------------------------------------------
  data->GetModuleLabel("mu2e::StepPointMCCollection",strh_module_label);
  data->GetDescription("mu2e::StepPointMCCollection",strh_description);

  art::Handle<mu2e::StepPointMCCollection>       strh_handle;
  const mu2e::StepPointMCCollection*             list_of_hits(0);

  static mu2e::GlobalConstantsHandle<mu2e::ParticleDataTable> pdt;
  mu2e::ParticleDataTable::maybe_ref info;
  
  if (strh_module_label[0] != 0) {
    if (strh_description[0] != 0) 
      AnEvent->getByLabel(strh_module_label, strh_description, strh_handle);
    if (strh_handle.isValid()) list_of_hits = strh_handle.product();
  }

  if (list_of_hits == NULL) {
    printf(" >>> ERROR in : StntupleInitMu2eVirtualDataBlock no list_of_hits. BAIL OUT\n");
    return -1;
  }

  //load time offset of this event
  fgTimeOffsets->updateMap(*AnEvent);

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  nhits = list_of_hits->size();

  const mu2e::StepPointMC* step;
  art::Ptr<mu2e::SimParticle> sim;

  TVdetHitData*           hit; 

  int                     vdIndex, pdg_id, gen_index;
  float time;
  float energyKin, energy;
  float mass;
  float mc_mom;
  float mc_momX, mc_momY, mc_momZ;
  float mc_posX, mc_posY, mc_posZ;

  for (int i=0; i<nhits; i++) {
    step  = &list_of_hits->at(i);

    sim = step->simParticle();
    if (!(sim->fromGenerator())) goto NEXT_VHIT;

    hit = data->NewHit();

    vdIndex   = step->volumeId();
    time      = fgTimeOffsets->timeWithOffsetsApplied(*step);

    pdg_id    = sim->pdgId();
    info      = pdt->particle(pdg_id);    
    
    mass      = info.ref().mass();
    energy    = sqrt(step->momentum().mag2() + std::pow(mass, 2));
    energyKin = energy - mass;

    gen_index = sim->genParticle()->generatorId().id();
          
    mc_mom    = step->momentum().mag();
    mc_momX   = step->momentum().x();
    mc_momY   = step->momentum().y();
    mc_momZ   = step->momentum().z();
//-----------------------------------------------------------------------------
// poor-man global-to-detector (tracker) coordinate system transformation
//-----------------------------------------------------------------------------
    mc_posX   = step->position().x()+3904.;
    mc_posY   = step->position().y();
    mc_posZ   = step->position().z()-10200;

    hit->Set(vdIndex, time, mass, energyKin, energy, 
	     pdg_id, gen_index, 
	     mc_mom, mc_momX, mc_momY, mc_momZ,
	     mc_posX, mc_posY, mc_posZ);
  NEXT_VHIT:;
  }

  data->f_RunNumber   = rn_number;
  data->f_EventNumber = ev_number;
  
  return 0;
}

