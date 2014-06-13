#ifndef murat_TStnClusterBlock
#define murat_TStnClusterBlock
//-----------------------------------------------------------------------------
//  definition of the cluster block for MU2E analysis
//  Author:    Pavel Murat (Fermilab)
//  Date:      March 07 2013
//-----------------------------------------------------------------------------
#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnCluster.hh"

class TStnClusterBlock: public TStnDataBlock {
  friend Int_t StntupleInitMu2eClusterBlockLinks(TStnDataBlock*, AbsEvent* , int);
public:
//----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  Int_t          fNClusters;
  TClonesArray*  fListOfClusters;
//----------------------------------------------------------------------------
//  functions
//----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStnClusterBlock();
  virtual ~TStnClusterBlock();

  TStnCluster* NewCluster() {
    TStnCluster* cl = new ((*fListOfClusters)[fNClusters]) TStnCluster(fNClusters);
    fNClusters++;
    return cl;
  }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Int_t           NClusters     () { return fNClusters;   }
  TClonesArray*   ListOfClusters() { return fListOfClusters; }

  TStnCluster*   Cluster(int I) {
    return (TStnCluster*) fListOfClusters->UncheckedAt(I); 
  }

//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* opt="") const;

  ClassDef(TStnClusterBlock,1)
};

#endif
