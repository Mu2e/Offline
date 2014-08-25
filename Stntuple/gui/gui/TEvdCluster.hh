///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdCluster_hh
#define TEvdCluster_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#include "Stntuple/base/TVisNode.hh"

#ifndef __CINT__
#include "RecoDataProducts/inc/CaloCluster.hh"
#else
namespace mu2e {
  class CaloCluster;
};
#endif


class TEvdCluster: public TObject {
public:
  
protected:

  TEllipse* fTrkEllipse;
  TEllipse* fCalEllipse;

  const mu2e::CaloCluster* fCluster;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdCluster() {}
  TEvdCluster(const mu2e::CaloCluster* fCluster); 

  virtual ~TEvdCluster();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::CaloCluster*  Cluster() { return fCluster; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  virtual void  Paint   (Option_t* Option = "");
  virtual void  PaintXY (Option_t* Option = "");
  virtual void  PaintCal(Option_t* Option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TEvdCluster,0)
};


#endif
