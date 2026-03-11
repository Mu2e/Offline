//
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"

#include "Offline/CalPatRec/inc/HlPrint.hh"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"


using namespace std;

namespace mu2e {

  HlPrint* HlPrint::_Instance(nullptr);

//-----------------------------------------------------------------------------
HlPrint::HlPrint(const fhicl::ParameterSet* PSet) {

  _event      = nullptr;
  //  _printUtils = new mu2e::TrkPrintUtils(PSet->get<fhicl::ParameterSet>("printUtils",fhicl::ParameterSet()));
}

//------------------------------------------------------------------------------
HlPrint* HlPrint::Instance(const fhicl::ParameterSet* PSet) {
  static HlPrint::Cleaner cleaner;

  if  (! _Instance) _Instance  = new HlPrint(PSet);
  return _Instance;
}


//______________________________________________________________________________
HlPrint::~HlPrint() {
  //  delete _printUtils;
}

//------------------------------------------------------------------------------
HlPrint::Cleaner::Cleaner() {
}


//------------------------------------------------------------------------------
  HlPrint::Cleaner::~Cleaner() {
    if (HlPrint::_Instance) {
      delete HlPrint::_Instance;
      HlPrint::_Instance = nullptr;
    }
  }


//-----------------------------------------------------------------------------
void HlPrint::printComboHit(const mu2e::ComboHit* Hit, const mu2e::StrawGasStep* Step,
                            const char* Opt,
                            int IHit, int Flags) {
  TString opt = Opt;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("#----------------------------------------------------------------------------------------------------");
    printf("--------------------------------------------------------------------------------------------\n");
    printf("#   I nsh   SID   Flags  Stn:Pln:Pnl:Str      X        Y        Z       Time     TCorr     eDep   End");
    printf("  DrTime  PrTime  TRes    WDist     WRes        PDG     PDG(M) GenID simID       p        pz\n");
    printf("#----------------------------------------------------------------------------------------------------");
    printf("--------------------------------------------------------------------------------------------\n");
  }

  if (opt == "banner") return;

  const mu2e::SimParticle * sim (0);

  int      pdg_id(-1), mother_pdg_id(-1), generator_id(-1), sim_id(-1);
  double   mc_mom(-1.);
  double   mc_mom_z(-1.);

  mu2e::GenId gen_id;

  if (Step) {
    art::Ptr<SimParticle> const& simptr = Step->simParticle();
    art::Ptr<SimParticle> mother        = simptr;

    while(mother->hasParent()) mother = mother->parent();

    sim           = mother.operator ->();

    pdg_id        = simptr->pdgId();
    mother_pdg_id = sim->pdgId();

    if (simptr->fromGenerator()) generator_id = simptr->genParticle()->generatorId().id();
    else                         generator_id = -1;

    sim_id        = simptr->id().asInt();
    mc_mom        = Step->momvec().mag();
    mc_mom_z      = Step->momvec().z();
  }

  if ((opt == "") || (opt.Index("data") >= 0)) {
    if (IHit  >= 0) printf("%5i " ,IHit);
    else            printf("      ");

    printf("%3i ",Hit->nStrawHits());

    printf("%5u",Hit->strawId().asUint16());

    if (Flags >= 0) printf(" %08x",Flags);
    else            printf("        ");

    printf(" %3i %3i %3i %3i",
           Hit->strawId().station(),
           Hit->strawId().plane(),
           Hit->strawId().panel(),
           Hit->strawId().straw());

    printf("  %8.3f %8.3f %9.3f %8.3f %8.3f %8.5f",
           Hit->pos().x(),Hit->pos().y(),Hit->pos().z(),
           Hit->time(),
           Hit->correctedTime(),
           Hit->energyDep());

    printf("  %3i %7.2f %7.2f %5.2f %8.3f %8.3f",
           (int) Hit->earlyEnd(),
           Hit->driftTime(),
           Hit->propTime(),
           Hit->transRes(),
           Hit->wireDist(),
           Hit->wireRes());

    printf(" %10i %10i %5i %5i %8.3f %8.3f\n",
           pdg_id,
           mother_pdg_id,
           generator_id,
           sim_id,
           mc_mom,
           mc_mom_z);
  }
}

//-----------------------------------------------------------------------------
// if 'ShfCollTag' is defined, print hit flags stored in the flag collection
// otherwise, print flags stored in the hit payload
//-----------------------------------------------------------------------------
void HlPrint::printComboHitCollection(const char* ChCollTag   ,
                                      const char* ShfCollTag  ,
                                      const char* SdmcCollTag ,
                                      double TMin, double TMax) {
//-----------------------------------------------------------------------------
// get combo hits
//-----------------------------------------------------------------------------
  art::Handle<mu2e::ComboHitCollection> chcH;
  const mu2e::ComboHitCollection*       chc(nullptr);

  _event->getByLabel(ChCollTag,chcH);

  if (chcH.isValid()) chc = chcH.product();
  else {
    printf("ERROR: cant find ComboHitCollection tag=%s, EXIT\n",ChCollTag);
    //    print_sh_colls();
    return;
  }

  art::Handle<mu2e::StrawHitFlagCollection> shfcH;
  const mu2e::StrawHitFlagCollection*       shfc(nullptr);
  if (*ShfCollTag != 0) {
    _event->getByLabel(ShfCollTag,shfcH);

    if (shfcH.isValid()) shfc = shfcH.product();
    else {
      printf("HlPrint::%s: ERROR: cant find StrawHitFlagCollection tag=%s, EXIT\n",__func__,ShfCollTag);
      // print_shf_colls();
      return;
    }
  }

  art::Handle<mu2e::StrawDigiMCCollection> sdmccH;
  const StrawDigiMCCollection*             sdmcc(nullptr);
  _event->getByLabel<mu2e::StrawDigiMCCollection>(SdmcCollTag,sdmccH);
  if (sdmccH.isValid())   sdmcc = sdmccH.product();
  else {
    printf("HlPrint::%s: ERROR: cant find StrawDigiMCCollection tag=%s, EXIT\n",__func__,SdmcCollTag);
    // print_sdmc_colls();
    return;
  }

  int nhits = chc->size();

  // const mu2e::ComboHit* hit0 = &chc->at(0);

  int banner_printed = 0;
  for (int i=0; i<nhits; i++) {
    const ComboHit* hit = &chc->at(i);
    int ind = hit->index(0);

    const mu2e::StrawDigiMC*  sdmc = &sdmcc->at(ind);
    const mu2e::StrawGasStep* step (nullptr);

    step = sdmc->earlyStrawGasStep().get();
//-----------------------------------------------------------------------------
// if SHF collection name is defined, print flags from that collection
// if not, print flags stored in the hit
// this is especially important for single-straw combo hits, as they are flagged
// only once - in MakeStrawHits module
//-----------------------------------------------------------------------------
    int flag(0);
    if (shfc) flag = *((int*) &shfc->at(i));
    else      flag = *((int*) &hit->flag());
    if (banner_printed == 0) {
      printComboHit(hit, step, "banner");
      banner_printed = 1;
    }
    if ((hit->time() >= TMin) && (hit->time() <= TMax)) {
      printComboHit(hit, step, "data", i, flag);
    }
  }

}

}
