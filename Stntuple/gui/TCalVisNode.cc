///////////////////////////////////////////////////////////////////////////////
// May 04 2013 P.Murat
// 
// in 'XY' mode draw calorimeter clusters as circles with different colors 
// in 'Cal' mode draw every detail...
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
#include "TObjArray.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// #include "Stntuple/geom/TTrajectory3D.hh"
// #include "Stntuple/geom/TTrajectoryPoint.hh"
// #include <TGeant/TG3Trap.hh>

// #include "murat/ana/TAnaDump.hh"

#include "Stntuple/gui/TCalVisNode.hh"
#include "Stntuple/gui/TEvdCalSection.hh"
#include "Stntuple/gui/TEvdCluster.hh"
#include "Stntuple/gui/TEvdCrystal.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "Stntuple/obj/TDisk.hh"

#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

ClassImp(TCalVisNode)

//_____________________________________________________________________________
  TCalVisNode::TCalVisNode(const char* Name, const mu2e::Disk* Disk, int   SectionID)
    :TVisNode(Name) 
{
  TEvdCrystal            *evd_cr;
  const mu2e::Crystal    *cr;

  fMinClusterEnergy = 5.;
  fMinCrystalEnergy = 0.;

  fSectionID         = SectionID;
  fListOfEvdCrystals = new TObjArray();
  fListOfEvdClusters = new TClonesArray("TEvdCluster",10);
  fEvdCalSection     = new TEvdCalSection(Disk, SectionID);

  fNCrystals = Disk->nCrystals();
					// assuming disks are the same
  fFirst     = SectionID*fNCrystals;

  TDisk* disk(0);
  
  for (int i=0; i<fNCrystals; i++) {
    cr     = &Disk->crystal(i);
    evd_cr = new TEvdCrystal(cr,disk);
    fListOfEvdCrystals->Add(evd_cr);
  }

  fTimePeak = NULL;
}

//_____________________________________________________________________________
TCalVisNode::~TCalVisNode() {
  fListOfEvdCrystals->Delete();
  delete fListOfEvdCrystals;

  delete fListOfEvdClusters;
  delete fEvdCalSection;
}


//-----------------------------------------------------------------------------
// 2014-08-11 : implement locally function removed from the mu2e::BaseCalorimeter
//-----------------------------------------------------------------------------
int TCalVisNode::LocalCrystalID(int CrystalID) {
  
  mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
  int nc, ns, id;

  ns = dc->nSection();
  id = CrystalID;

  for (int i=0; i<ns; ++i) {
    nc = dc->section(i).nCrystals();
    if (id < nc) break;
    id -= nc;
  }
  return id;
}

//-----------------------------------------------------------------------------
//  each CAL vis node corresponds to one disk with a given fSectionID
//-----------------------------------------------------------------------------
int TCalVisNode::InitEvent() {

  art::ServiceHandle<mu2e::GeometryService> geom;

  mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;

  //  static int first_call(1);

  const mu2e::CaloCluster*     cl;
  const mu2e::CaloCrystalHit*  hit;

  TEvdCluster*              evd_cl;
  TEvdCrystal*              evd_cr;

  int      ncl, loc, nhits, id;
  double   energy;
//-----------------------------------------------------------------------------
// individual crystals, assume that both disks have the same size
//-----------------------------------------------------------------------------
  Clear();
  
  if (*fListOfCrystalHits != NULL) {
    nhits = (*fListOfCrystalHits)->size();
    for (int i=0; i<nhits; i++) {
      hit = &(*fListOfCrystalHits)->at(i);
					// in short, the crystal number
      id  = hit->id()-fFirst;
      if ((id >= 0) && (id < fNCrystals)) {
//-----------------------------------------------------------------------------
// hit on a given disk
//-----------------------------------------------------------------------------
	if (hit->energyDep() > fMinCrystalEnergy) {
	  evd_cr = EvdCrystal(id);
	  
	  evd_cr->AddHit(hit);
	  evd_cr->SetFillColor(kYellow-7); 
	  evd_cr->SetFillStyle(1024);
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// clusters
//-----------------------------------------------------------------------------
  if (*fListOfClusters != NULL) {

    ncl = (*fListOfClusters)->size();

    for (int i=0; i<ncl; i++) {
      cl = &(*fListOfClusters)->at(i);
      if (cl->sectionId() == fSectionID) {
	evd_cl = NewEvdCluster(cl);
//-----------------------------------------------------------------------------
// set colors of the crystals
//-----------------------------------------------------------------------------
	const mu2e::CaloCluster::CaloCrystalHitPtrVector caloClusterHits = cl->caloCrystalHitsPtrVector();
	int nh = caloClusterHits.size();

	for (int i=0; i<nh; i++) {
	  const mu2e::CaloCrystalHit* hit = &(*caloClusterHits.at(i));
	  int id = hit->id();

	  loc = LocalCrystalID(id);
//-----------------------------------------------------------------------------
// find a crystal with a given ID, display it in red
//-----------------------------------------------------------------------------
	  evd_cr = (TEvdCrystal*) fListOfEvdCrystals->At(loc);

	  evd_cl->AddCrystal(evd_cr);
//-----------------------------------------------------------------------------
// displayed color of the crystal is define by the max hit energy
//-----------------------------------------------------------------------------
	  energy = hit->energyDep();

	  if (energy > fMinCrystalEnergy) {
	    if      (energy > 100.) evd_cr->SetFillColor(kRed+ 2);
	    else if (energy >  10.) evd_cr->SetFillColor(kRed   ); 
	    else if (energy >   1.) evd_cr->SetFillColor(kRed- 9); 
	    else                    evd_cr->SetFillColor(kRed-10); 
	    
	    evd_cr->SetFillStyle(1024);
	  }
	}
      }
    }
  }

  return 0;
}



//_____________________________________________________________________________
void TCalVisNode::Paint(Option_t* Option) {
  // paints one disk (.. or vane, in the past), i.e. section

  int   iv;
				// parse option list
  const char* view = TVisManager::Instance()->GetCurrentView();

  if      (strstr(view,"trkxy" ) != 0) {
    PaintXY (Option);
  }
  if      (strstr(view,"cal"   ) != 0) {
    sscanf(view,"cal,%i",&iv);
    if (iv == fSectionID) {
      PaintCal(Option);
    }
  }
  else {
    // what is the default?
    //    Warning("Paint",Form("Unknown option %s",option));
  }

  gPad->Modified();
}


//-----------------------------------------------------------------------------
// I have 2 vis nodes, one for each disk
//-----------------------------------------------------------------------------
void TCalVisNode::PaintXY(Option_t* Option) {
  // draw calorimeter
 
  fEvdCalSection->PaintXY(Option);
//-----------------------------------------------------------------------------
// now can draw clusters
//-----------------------------------------------------------------------------
  TEvdCluster* cl;
  TString      opt;
  int          display_cluster;
  double       time;
   
  TStnVisManager* vm = TStnVisManager::Instance();

  int ipeak = vm->TimePeak();

  if (ipeak >= 0) {
    if ((*fCalTimePeakColl) != NULL) {
      int ntp = (*fCalTimePeakColl)->size();
      if (ipeak < ntp) fTimePeak = &(*fCalTimePeakColl)->at(ipeak);
      else             fTimePeak = NULL;
    }
  }

  int ncl = fListOfEvdClusters->GetEntries();
  if (ncl > 0) {
    //    TAnaDump::Instance()->printCaloCluster(0,"banner");
    for (int i=0; i<ncl; i++) {
      cl = EvdCluster(i);
      //      TAnaDump::Instance()->printCaloCluster(cl,"data");

      time = cl->Cluster()->time();
      display_cluster = 1;

      if (fTimePeak != NULL) {
	if ((time < fTimePeak->TMin()) || (time > fTimePeak->TMax())) display_cluster = 0;
      }

      if (display_cluster) {
//-----------------------------------------------------------------------------
// display only clusters with E > 5 MeV
//-----------------------------------------------------------------------------
	if (cl->Cluster()->energyDep() > fMinClusterEnergy) {
	  cl->Paint(Option);
	}
      }
    }
  }
}

//_____________________________________________________________________________
void TCalVisNode::PaintCal(Option_t* Option) {
  // draw calorimeter

  fEvdCalSection->PaintCal(Option);

  art::ServiceHandle<mu2e::GeometryService> geom;
  mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;

//-----------------------------------------------------------------------------
// draw clusters ... with the size proportional to what?  
//-----------------------------------------------------------------------------
  TEvdCluster* cl;

  int ncl = fListOfEvdClusters->GetEntries();

  for (int i=0; i<ncl; i++) {
    cl    = EvdCluster(i);
//-----------------------------------------------------------------------------
// draw crystals
//-----------------------------------------------------------------------------
    TEvdCrystal*             cr;

    int ncr  = NCrystals();

    for (int i=0; i<ncr; i++) {
      cr = EvdCrystal(i);
      cr->SetLineWidth(0.1);
      cr->PaintCal(Option);
    }
//-----------------------------------------------------------------------------
// display only clusters with E > 5 MeV
// crystals included into these clusters have thick red border
//-----------------------------------------------------------------------------
    if (cl->Cluster()->energyDep() > fMinClusterEnergy) {

      int ncc = cl->NCrystals();
      for (int i=0; i<ncc; i++) {
	cr = cl->Crystal(i);
	cr->SetLineColor(kRed);
	cr->SetLineWidth(2);
	cr->PaintCal(Option);
      }
      cl->PaintCal(Option);
    }
  }
}


//_____________________________________________________________________________
void TCalVisNode::PaintRZ(Option_t* option) {
  // draw calorimeter
}

//_____________________________________________________________________________
Int_t TCalVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
  int              min_dist(9999), iv;
  static TVector3  global;

  TVisManager* vm = TVisManager::Instance();
  const char* view = vm->GetCurrentView();

  if      (strstr(view,"trkxy" ) != 0) return min_dist;
  
  sscanf(view,"cal,%i",&iv);
  if (iv != fSectionID)                return min_dist;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  //  printf("px,py,X,Y = %5i %5i %10.3f %10.3f\n",px,py,global.X(),global.Y());

  min_dist = DistancetoPrimitiveXY(px,py);

  return min_dist;
}

//_____________________________________________________________________________
Int_t TCalVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t min_dist = 9999;

  static TVector3 global;
  static TVector3 local;
  int             px0, py0;
  double          dist;

  fClosestObject = NULL;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  int          ncr;
  TEvdCrystal  *cr;

  ncr = fListOfEvdCrystals->GetEntries();
//   x   = global.X();
//   y   = global.Y();

  for (int icr=0; icr<ncr; icr++) {
    cr   = EvdCrystal(icr);
    
    px0 = gPad->XtoAbsPixel(cr->X0());
    py0 = gPad->YtoAbsPixel(cr->Y0());
    dist = sqrt((px0-px)*(px0-px) + (py0-py)*(py0-py));
    if (dist < min_dist) {
      min_dist       = dist;
      fClosestObject = cr;
    }
  }

  return min_dist;
}

//_____________________________________________________________________________
Int_t TCalVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}


//-----------------------------------------------------------------------------
void TCalVisNode::Clear(Option_t* Opt) {
  int          ncr;
  TEvdCrystal  *evd_cr;

  ncr = fListOfEvdCrystals->GetEntries();
  for (int icr=0; icr<ncr; icr++) {
    evd_cr = EvdCrystal(icr);
    evd_cr->Clear();
  }

  fListOfEvdClusters->Clear();
  fNClusters = 0;
}

//-----------------------------------------------------------------------------
void TCalVisNode::Print(Option_t* Opt) const {
  printf(" >>> TCalVisNode::Print is not implemented yet\n");
}

