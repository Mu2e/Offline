///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

#include "Stntuple/obj/TGenpBlock.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"

//_____________________________________________________________________________
int StntupleInitMu2eGenpBlock(TStnDataBlock* block, AbsEvent* AnEvent, int mode) 
{
  // fill generator particle data block (GENP - GENerator Particles - 
  // is the name from Run I)

  std::vector<art::Handle<mu2e::GenParticleCollection>> list_of_gp;

  const mu2e::GenParticleCollection*     coll(0);
  const mu2e::GenParticle*               gp(0);
  //  const art::Provenance*                 prov;
  const art::Handle<mu2e::GenParticleCollection>* handle;

  art::Selector  selector(art::ProductInstanceNameSelector(""));

  double px, py, pz, mass, energy;
  int    pdg_code;
    //  TLorentzVector v;

  TGenpBlock* genp_block = (TGenpBlock*) block;
  genp_block->Clear();
//-----------------------------------------------------------------------------
// initialization from HEPG, loop over HEPG particles and fill the list
// loop over the existing HEPG banks, add non-initialized interaction
//-----------------------------------------------------------------------------
  AnEvent->getMany(selector, list_of_gp);

  TDatabasePDG* pdg_db = TDatabasePDG::Instance();
  TParticlePDG* part;

  for (std::vector<art::Handle<mu2e::GenParticleCollection>> ::const_iterator it = list_of_gp.begin();
       it != list_of_gp.end(); it++) {
    handle = it.operator -> ();

    if (handle->isValid()) {
      coll = handle->product();
      //      prov = handle->provenance();

//       printf("moduleLabel = %-20s, producedClassname = %-30s, productInstanceName = %-20s\n",
// 	     prov->moduleLabel().data(),
// 	     prov->producedClassName().data(),
// 	     prov->productInstanceName().data());

      for (std::vector<mu2e::GenParticle>::const_iterator ip = coll->begin();
	   ip != coll->end(); ip++) {
	gp       = ip.operator -> ();
	pdg_code = (int) gp->pdgId();
	part     = pdg_db->GetParticle(pdg_code);

	px     = gp->momentum().x();
	py     = gp->momentum().y();
	pz     = gp->momentum().z();
	mass   = part->Mass();
	energy = sqrt(px*px+py*py+pz*pz+mass*mass);

	genp_block->NewParticle(pdg_code,
				(int) gp->generatorId().id(),
				-1, -1, -1, -1,
				px, py, pz, energy,
				gp->position().x(),
				gp->position().y(),
				gp->position().z(),
				gp->time(),
				gp->properTime());
      }
    }
    else {
      printf(">>> ERROR in StntupleInitMu2eGenpBlock: failed to locate collection");
      printf(". BAIL OUT. \n");
      return -1;
    }
  }

  return 0;
}


