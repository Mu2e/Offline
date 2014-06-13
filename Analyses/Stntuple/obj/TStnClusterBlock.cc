#include <iostream>
#include <iomanip>

#include "obj/TStnClusterBlock.hh"
#include "obj/TStnCluster.hh"

ClassImp(TStnClusterBlock)
//______________________________________________________________________________
void TStnClusterBlock::Streamer(TBuffer &R__b) {
  // Stream an object of class TStnClusterBlock.

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    R__b >> fNClusters;
    fListOfClusters->Streamer(R__b);
    for (int i=0; i<fNClusters; i++) {
      Cluster(i)->SetNumber(i);
    }
  } 
  else {
    R__b.WriteVersion(TStnClusterBlock::IsA());
    R__b << fNClusters;
    fListOfClusters->Streamer(R__b);
  }
}

//_____________________________________________________________________________
TStnClusterBlock::TStnClusterBlock() {
  fNClusters   = 0;
  fListOfClusters = new TClonesArray("TStnCluster",100);
  fListOfClusters->BypassStreamer(kFALSE);
  fCollName  = "default";
}


//_____________________________________________________________________________
TStnClusterBlock::~TStnClusterBlock() {
  fListOfClusters->Delete();
  delete fListOfClusters;
}


//_____________________________________________________________________________
void TStnClusterBlock::Clear(Option_t* opt) {
  fNClusters = 0;
  fListOfClusters->Clear(opt);
}

//------------------------------------------------------------------------------
void TStnClusterBlock::Print(Option_t* opt) const {

  int banner_printed = 0;
  for (int i=0; i<fNClusters; i++) {
    TStnCluster* t = ((TStnClusterBlock*) this)->Cluster(i);
    if (! banner_printed) {
      t->Print("banner");
      banner_printed = 1;
    }
    t->Print("data");
  }
}
