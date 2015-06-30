///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include "Stntuple/obj/TSimpBlock.hh"


ClassImp(TSimpBlock)

//______________________________________________________________________________
void TSimpBlock::Streamer(TBuffer &R__b) 
{
   // Stream an object of class TSimpBlock as compact, as possible

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    R__b >> fNParticles;
    fListOfParticles->Streamer(R__b);
    for (int i=0; i<fNParticles; i++) {
      Particle(i)->SetUniqueID(i);
    }
  } 
  else {
    R__b.WriteVersion(TSimpBlock::IsA());
    R__b << fNParticles;
    fListOfParticles->Streamer(R__b);
  }
}

//_____________________________________________________________________________
TSimpBlock::TSimpBlock() {
  fNParticles      = 0;
  fListOfParticles = new TClonesArray("TSimParticle",10);
  fListOfParticles->BypassStreamer(kFALSE);
}


//_____________________________________________________________________________
TSimpBlock::~TSimpBlock() {
  fListOfParticles->Delete();
  delete fListOfParticles;
}


//_____________________________________________________________________________
void TSimpBlock::Clear(const char* opt) {
  fNParticles    = 0;
					// don't modify cut values at run time
  fListOfParticles->Clear(opt);
}

//_____________________________________________________________________________
TSimParticle* 
TSimpBlock::NewParticle(Int_t ID, Int_t ParentID, Int_t PdgCode, 
			int CreationCode, int TerminationCode,
			int StartVolumeIndex, int EndVolumeIndex,
			Float_t px, Float_t py, Float_t pz, Float_t e,
			Float_t vx, Float_t vy, Float_t vz, Float_t t)
{
  // add new particle to the block. Block is filled sequentially, so 
  // assume that this is the last particle (not necessarily!) and it is 
  // added to the last primary interaction

  TSimParticle* p;

  p = new ((*fListOfParticles)[fNParticles]) TSimParticle();

  p->Init(ID,ParentID, PdgCode,CreationCode,TerminationCode,
	  StartVolumeIndex,EndVolumeIndex,
	  px,py,pz,e,vx,vy,vz,t);

  fNParticles += 1;

  return p;
}


//_____________________________________________________________________________
void TSimpBlock::Print(const char* Opt) const {
  // opt: /c : comment lines only, useful for printing the hard interaction only

  int banner_printed(0), i1, i2 /*, comments_only(0)*/;

  //  if (strstr(Opt,"/c") != 0) comments_only = 1;

  TSimpBlock* simp = (TSimpBlock*) this;

  i1 = 0;
  i2 = simp->NParticles();

  for (int i=i1; i<i2; i++) {
    TSimParticle* p = simp->Particle(i);
    if (! banner_printed) {
      p->Print("banner");
      banner_printed = 1;
    }
    p->Print("data");
  }
}


