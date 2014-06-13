///////////////////////////////////////////////////////////////////////////////
// 2014-01-26 P.Murat
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
//_____________________________________________________________________________
Int_t StntupleInitMu2eStrawDataBlock(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) 
{
  // initialize COT data block with the `event' data
  // return -1, if bank doesn't exist, 0, if everything is OK

  int ev_number, rn_number, nhits;

  static char   strh_module_label[100], strh_description[100];

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

  TStrawDataBlock* data = (TStrawDataBlock*) Block;
  data->Clear();
//-----------------------------------------------------------------------------
//  straw hit information
//-----------------------------------------------------------------------------
  data->GetModuleLabel("mu2e::StrawHitCollection",strh_module_label);
  data->GetDescription("mu2e::StrawHitCollection",strh_description );

  art::Handle<mu2e::StrawHitCollection> strh_handle;
  const mu2e::StrawHitCollection*             list_of_hits(0);

  if (strh_module_label[0] != 0) {
    if (strh_description[0] == 0) AnEvent->getByLabel(strh_module_label,strh_handle);
    else                          AnEvent->getByLabel(strh_module_label,strh_description,strh_handle);
    list_of_hits = /*(mu2e::StrawHitCollection*) */ strh_handle.product();
  }

  if (list_of_hits == NULL) {
    printf(" >>> ERROR in StntupleInitMu2eCalDataBlock: no list_of_hits. BAIL OUT\n");
    return -1;
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  nhits = list_of_hits->size();

  const mu2e::StrawHit* sh; 
  TStrawHitData*        hit; 

  for (int i=0; i<nhits; i++) {
    sh  = &list_of_hits->at(i);
    hit = data->NewHit();
    hit->Set(sh->strawIndex().asInt(),sh->time(),sh->dt(),sh->energyDep());
  }

  data->f_RunNumber = rn_number;
  data->f_EventNumber = ev_number;
  
  return 0;
}

