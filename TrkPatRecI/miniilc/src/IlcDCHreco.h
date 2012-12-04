#ifndef ILCDCHRECO_H
#define ILCDCHRECO_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                   DCH reconstruction name space
//
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------
#include <Rtypes.h>

//namespace IlcDCHreco {    
   const Int_t kMaxClusterPerLayer=2000;
   const Int_t kRowsToSkip=10;
   const Int_t kMaxLayer=40;
   const Int_t kMaxInLayer=5;

   const Double_t kMaxCHI2cluster=16.;
   const Double_t kMaxCHI2track=3.;
   const Double_t kMaxROAD=30.;

namespace ClusterIndex{
  enum Index_t{
      Sign   =0x80000000,
      Sector =0xff000000,
      Layer  =0x00fff000,
      NotAccept =0x00000800,
      Cluster=0x000007ff,
      SignOffset   =31,
      SectorOffset =24,
      LayerOffset  =12,
      NotAcceptOffset =11,
      ClusterOffset=0
  };
}
//}

//using namespace IlcDCHreco;

#endif
