///////////////////////////////////////////////////////////////////////////////
// May 04 2013 P.Murat
// 
///////////////////////////////////////////////////////////////////////////////
#include "TVirtualX.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArc.h"
#include "TArrow.h"
#include "TMath.h"
#include "TBox.h"
#include "TEllipse.h"
#include "TObjArray.h"

// #include "art/Framework/Principal/Handle.h"

// #include "GeometryService/inc/GeometryService.hh"
// #include "GeometryService/inc/GeomHandle.hh"
// #include "CalorimeterGeom/inc/VaneCalorimeter.hh"

// #include "CalorimeterGeom/inc/Crystal.hh"
// #include "CalorimeterGeom/inc/Disk.hh"
// #include "CalorimeterGeom/inc/DiskCalorimeter.hh"
// #include "CalorimeterGeom/inc/Calorimeter.hh"

#include "Stntuple/obj/TDisk.hh"
#include "Stntuple/obj/TDiskCalorimeter.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStnCrystal.hh"

ClassImp(TDiskCalorimeter)

//_____________________________________________________________________________
TDiskCalorimeter::TDiskCalorimeter(): TObject() {
  fInitialized = 0;
}

//_____________________________________________________________________________
TDiskCalorimeter::TDiskCalorimeter(GeomData_t* Geom): TObject() {
  fInitialized = 0;
  Init(Geom);
}

//-----------------------------------------------------------------------------
TDiskCalorimeter::~TDiskCalorimeter() {
  for (int i=0; i<fNDisks; i++) {
    delete fDisk[i];
  }
}

//-----------------------------------------------------------------------------
int TDiskCalorimeter::Init(GeomData_t* Geom) {
  int channel_offset(0);

  if (fInitialized != 0) {
    printf(">>> TDiskCalorimeter::Init ERROR: an attempt to reinitialize, BAIL OUT\n");
    return -1;
  }

  fNDisks = Geom->fNDisks;

  double dead_space = Geom->fWrapperThickness+Geom->fShellThickness;

  for (int i=0; i<fNDisks; i++) {
    fDisk[i] = new TDisk(i,
			 Geom->fRMin[i],
			 Geom->fRMax[i],
			 Geom->fZ0  [i],
			 Geom->fHexSize,
			 dead_space,
			 Geom->fMinFraction);

    fDisk[i]->SetChannelOffset(channel_offset);
    channel_offset += fDisk[i]->NCrystals();
  }

  fInitialized = 1;

  return 0;
}

//-----------------------------------------------------------------------------
// so far assume that 'HitID' is just a sequential channel number
//-----------------------------------------------------------------------------
int TDiskCalorimeter::DiskNumber(int HitID) {
  int      dn(-1), first, last;
  TDisk*   disk;

  for (int i=0; i<fNDisks; i++) {
    disk = fDisk[i];
    first = disk->ChannelOffset();
    last  = first+disk->NCrystals()-1;
    if ((HitID >= first) && (HitID <= last)) {
      dn = i;
      break;
    }
  }
  return dn;
}
//-----------------------------------------------------------------------------
// so far assume that 'HitID' is just a sequential channel number
//-----------------------------------------------------------------------------
double TDiskCalorimeter::CrystalRadius(int HitID) {
  int     dn, offset;
  double  r;

  dn = DiskNumber(HitID);

  TDisk* disk = Disk(dn);

  offset = HitID-disk->ChannelOffset();

  r = disk->Crystal(offset)->Radius();
 
  return r;
}

//-----------------------------------------------------------------------------
int TDiskCalorimeter::InitEvent(TCalDataBlock* CalDataBlock) {

  int               nhits, hit_id, dn, rc(0), loc;
  TStnCrystal*      crystal;
  TCalHitData*      hit;
  TDisk*            disk;

  Clear();

  nhits = CalDataBlock->NHits();
  
  for (int i=0; i<nhits; i++) {
    hit     = CalDataBlock->CalHitData(i);
    hit_id  = hit->ID();

    dn      = DiskNumber(hit_id);
    disk    = Disk(dn);
    loc     = hit_id-disk->ChannelOffset();
    crystal = disk->Crystal(loc);

    if (crystal != NULL) {
      crystal->AddHit(hit);
    }
    else {
      printf(">>> ERROR in TDiskCalorimeter::InitEvent : hit_id=%5i not assigned\n",hit_id);
      rc = -1;
    }
  }
  return rc;
}

//-----------------------------------------------------------------------------
void TDiskCalorimeter::Clear(Option_t* Opt) {
  for (int i=0; i<fNDisks; i++) {
    fDisk[i]->Clear();
  }
}

//-----------------------------------------------------------------------------
void TDiskCalorimeter::Print(Option_t* Opt) const {
  printf(">>> ERROR: TDiskCalorimeter::Print not implemented yet aa\n");
  
}


