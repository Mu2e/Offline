//-----------------------------------------------------------------------------
//  Apr 2013 P.Murat: initialization of the MU2E STNTUPLE cluster block
//
//-----------------------------------------------------------------------------
#include <cstdio>
#include "TROOT.h"
#include "TFolder.h"
#include "TLorentzVector.h"

#include "Stntuple/obj/TStnDataBlock.hh"

#include "Stntuple/obj/TStnCluster.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "KalmanTests/inc/KalRepPtrCollection.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"

// #include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"
// #include "TrackCaloMatching/inc/TrackClusterLink.hh"
// #include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

//-----------------------------------------------------------------------------
// assume that the collection name is set, so we could grab it from the event
//-----------------------------------------------------------------------------
int  StntupleInitMu2eClusterBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode) {

  const char*               oname = {"MuratInitClusterBlock"};
  
//   int                           station, ntrk;
//   KalRep                        *krep;  
//  double                        h1_fltlen, hn_fltlen, entlen, fitmom_err;
//   TStnTrack*                    track;
//   const mu2e::StepPointMC*      step;
  mu2e::CaloClusterCollection*  list_of_clusters;

  const double kMinECrystal = 0.1; // count crystals above 100 KeV

  static char                calo_module_label[100], calo_description[100]; 
  static char                trcl_module_label[100], trcl_description[100];

  TStnClusterBlock*         cb = (TStnClusterBlock*) Block;
  TStnCluster*              cluster;

  //  static int  first_entry(1);

  cb->Clear();

  //
  // "makeCaloCluster" would be the process name, "AlgoCLOSESTSeededByENERGY" - the description,
  // 
  
  cb->GetModuleLabel("mu2e::CaloClusterCollection",calo_module_label);
  cb->GetDescription("mu2e::CaloClusterCollection",calo_description);

  cb->GetModuleLabel("mu2e::TrackClusterLink",trcl_module_label);
  cb->GetDescription("mu2e::TrackClusterLink",trcl_description );

  art::Handle<mu2e::CaloClusterCollection> calo_cluster_handle;
  if (calo_description[0] == 0) Evt->getByLabel(calo_module_label,calo_cluster_handle);
  else                          Evt->getByLabel(calo_module_label,calo_description,calo_cluster_handle);
  list_of_clusters = (mu2e::CaloClusterCollection*) &(*calo_cluster_handle);


//   art::Handle<mu2e::TrackClusterLink>  trk_cal_map;
//   if (trcl_module_label[0] != 0) {
//     if (trcl_description[0] != 0) Evt->getByLabel(trcl_module_label,trcl_description,trk_cal_map);
//     else                          Evt->getByLabel(trcl_module_label,trk_cal_map);
//   }

  art::ServiceHandle<mu2e::GeometryService> geom;

  const mu2e::Calorimeter* cal;

  if      (geom->hasElement<mu2e::Calorimeter>() ) {
    mu2e::GeomHandle<mu2e::Calorimeter> cc;
    cal = cc.operator->();
  }
  else if (geom->hasElement<mu2e::VaneCalorimeter>() ) {
    mu2e::GeomHandle<mu2e::VaneCalorimeter> vc;
    cal = vc.operator->();
  }
  else if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
    mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
    cal = dc.operator->();
  }
//-----------------------------------------------------------------------------
// tracks are supposed to be already initialized
//-----------------------------------------------------------------------------
  const mu2e::CaloCluster       *cl;
  const mu2e::CaloCrystalHit    *hit;
  int                           id, ncl;

  double                        sume, sume2, sumy, sumz, sumy2, sumz2, qn;
  double                        e, e1, e2, emean, e2mean;
  double                        ymean, zmean, y2mean, z2mean;
  
  ncl = list_of_clusters->size();
  for (int i=0; i<ncl; i++) {
    cluster               = cb->NewCluster();
    cl                    = &list_of_clusters->at(i);
    cluster->fCaloCluster = cl;
    cluster->fDiskID      = cl->vaneId();
    cluster->fEnergy      = cl->energyDep();
    cluster->fTime        = cl->time();
    
    const mu2e::CaloCrystalHitPtrVector list_of_crystals = cluster->fCaloCluster->caloCrystalHitsPtrVector();

    int nh = list_of_crystals.size();
    Hep3Vector pos;
//-----------------------------------------------------------------------------
// print individual crystals in local vane coordinate system
// Y and Z 
//-----------------------------------------------------------------------------
    qn     = 0;
    sume   = 0;
    sume2  = 0;
    sumy   = 0;
    sumz   = 0;
    sumy2  = 0;
    sumz2  = 0;
    
    e2     = 0;
    
    for (int ih=0; ih<nh; ih++) {
      hit = &(*list_of_crystals.at(ih));
      e   = hit->energyDep();
      id  = hit->id();
      pos = cal->crystalOriginInSection(id);
      
      if (e > kMinECrystal) {
	qn    += 1.;
	sume  += e;
	sume2 += e*e;
	sumy  += e*pos.y();
	sumz  += e*pos.z();
	sumy2 += e*pos.y()*pos.y();
	sumz2 += e*pos.z()*pos.z();
	
	if (ih<2) {
	  e2 += e;
	  if (ih == 0) e1 = e;
	}
      }
    }
    
    emean  = sume/(qn+1.e-12);
    e2mean = sume2/(qn+1.e-12);
    
    ymean  = sumy /(sume+1.e-12);
    zmean  = sumz /(sume+1.e-12);
    y2mean = sumy2/(sume+1.e-12);
    z2mean = sumz2/(sume+1.e-12);

    cluster->fX         = cl->cog3Vector().x();
    cluster->fY         = cl->cog3Vector().y();
    cluster->fZ         = cl->cog3Vector().z();
    cluster->fIx1       = cl->cogRow();
    cluster->fIx2       = cl->cogColumn();

    cluster->fNCrystals = nh;
    cluster->fNCr1      = qn;
    cluster->fYMean     = ymean;
    cluster->fZMean     = zmean;
    cluster->fSigY      = sqrt(y2mean-ymean*ymean);
    cluster->fSigZ      = sqrt(z2mean-zmean*zmean);
    cluster->fSigR      = sqrt(y2mean+z2mean-ymean*ymean-zmean*zmean);
//-----------------------------------------------------------------------------
// make sure that nothing goes into an overflow
//_____________________________________________________________________________
    cluster->fFrE1      = e1/(sume+1.e-5);
    cluster->fFrE2      = e2/(sume+1.e-5);
    cluster->fSigE1     = sqrt(qn*(e2mean-emean*emean))/(sume+1.e-12);
    cluster->fSigE2     = sqrt(qn*(e2mean-emean*emean))/(emean+1.e-12);
    
//     unsigned int nm = (*trk_cal_map).size();
//     for(size_t im=0; i<nm; im++) {
//       //	KalRepPtr const& trkPtr = fTrkCalMap->at(i).first->trk();
//       //	const KalRep *  const &trk = *trkPtr;
      
//       cl = &(*(*trk_cal_map).at(im).second);
      
//       if (cl == cluster->fCaloCluster) { 
// 	cluster->fClosestTrack = fTrackBlock->Track(im);
// 	break;
//       }
//     }
  }
  return 0;
}

//_____________________________________________________________________________
Int_t StntupleInitMu2eClusterBlockLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) 
{
  // Mu2e version, do nothing

  Int_t  ev_number, rn_number;

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (! Block->Initialized(ev_number,rn_number)) return -1;

					// do not do initialize links 2nd time

  if (Block->LinksInitialized()) return 0;

  TStnClusterBlock* header = (TStnClusterBlock*) Block;
  //  TStnEvent* ev   = header->GetEvent();
//-----------------------------------------------------------------------------
// mark links as initialized
//-----------------------------------------------------------------------------
  header->fLinksInitialized = 1;

  return 0;
}

