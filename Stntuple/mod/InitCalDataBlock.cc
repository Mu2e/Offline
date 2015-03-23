///////////////////////////////////////////////////////////////////////////////
// 2014-01-26 P.Murat Mu2e version of TCalDataBlock initialization
///////////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include "TROOT.h"
#include "TFolder.h"
#include "TLorentzVector.h"

#include "Stntuple/obj/TCalDataBlock.hh"

#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
//-----------------------------------------------------------------------------
Int_t StntupleInitMu2eCalDataBlock(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) 
{
  // initialize CAL data block with the `event' data

  static char   calo_module_label[100], calo_description[100];
  int           ev_number, rn_number, nhits;

  mu2e::CaloCrystalHitCollection* list_of_hits(0);

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

  TCalDataBlock* data = (TCalDataBlock*) Block;
  data->Clear();

  data->GetModuleLabel("mu2e::CaloCrystalHitCollection",calo_module_label);
  data->GetDescription("mu2e::CaloCrystalHitCollection",calo_description );

       // Get handles to calorimeter crystal hits

  art::Handle<mu2e::CaloCrystalHitCollection> calo_hits_handle;

  if (calo_module_label[0] != 0) {
    if (calo_description[0] == 0) {
      AnEvent->getByLabel(calo_module_label,calo_hits_handle);
    }
    else {
      AnEvent->getByLabel(calo_module_label,
			  calo_description,
			  calo_hits_handle);
    }

    list_of_hits  = (mu2e::CaloCrystalHitCollection*) calo_hits_handle.product();
  }

  if (list_of_hits == NULL) {
    printf(" >>> ERROR in StntupleInitMu2eCalDataBlock: no list_of_hits. BAIL OUT\n");
    return -1;
  }

  nhits = list_of_hits->size();

  mu2e::CaloCrystalHit* calo_hit;
  TCalHitData*          hit;

  // reminder: data->fNHits is set to 0 by TCalDataBlock::Clear(), should be this way

  for (int i=0; i<nhits; i++) {
    calo_hit = &list_of_hits->at(i);
    hit      = data->NewCalHitData();

    hit->Set(calo_hit->id(),
	     calo_hit->numberOfROIdsUsed(),
	     calo_hit->time(),
	     calo_hit->energyDep());
  }
//-----------------------------------------------------------------------------
// store geometry data (this is, obviously, a kludge)
//-----------------------------------------------------------------------------
  art::ServiceHandle<mu2e::GeometryService> geom;
  mu2e::GeomHandle<mu2e::DiskCalorimeter>   dc;

  const mu2e::DiskCalorimeter*              cal;
  const mu2e::Disk*                         disk;

  cal = dc.get();

  data->fNDisks = cal->nDisk();
  for (int i=0; i<data->fNDisks; i++) {
    disk = &cal->disk(i);
    data->fRMin[i] = disk->innerRadius();
    data->fRMax[i] = disk->outerRadius();
    data->fZ0  [i] = disk->origin().z();
  }
  data->fCrystalSize = cal->caloGeomInfo().crystalHalfTrans();

				        // also a dummy line
  data->fMinFraction      = 1.0;
  data->fWrapperThickness = cal->caloGeomInfo().wrapperThickness();
  data->fShellThickness   = cal->caloGeomInfo().caseThickness  ();

					// on return set event and run numbers
					// to mark block as initialized
  data->f_RunNumber   = rn_number;
  data->f_EventNumber = ev_number;

  return 0;
}







