///////////////////////////////////////////////////////////////////////////////
/*
.L help.C
.L conv_ana.C
help()
conv_ana("stntuple/dev_240","tmet08-mj07","conv/debug/grl=etf,141544,168889",1)
x(151978,1201384)
TStntuple::SwimToCes(m_cnv->fExtrapolator,m_cnv->fTrackBlock->Track(0),1);

 */

// #include "global_vars.h"

// TStnTrackBlock*      trkb;
// TStnTauBlock*        taub;
// TStnPi0Block*        pib;
// TStnJetBlock*        jb;
// TStnConversion*      cnv;
// TStnElectronBlock*   eb;
// TStnMuonBlock*       mb;
// TStnConversionBlock* cnvb;

// //_____________________________________________________________________________
// void next_event() {
//   g.x->Continue(1);
// }

// //_____________________________________________________________________________
// TTree* get_stntuple(const char* Filename) {
//   TFile* f = TFile::Open(Filename);
//   TTree* t = (TTree*) f->Get("STNTUPLE");
//   return t;
// }

//_____________________________________________________________________________
void print_StrawHitCollection() {
  printf("\n");
  TAnaDump* d = TAnaDump::Instance();
  d->printStrawHitCollection("makeSH");
}

//_____________________________________________________________________________
void print_StepPointMCCollection() {
  printf("\n");
  TAnaDump* d = TAnaDump::Instance();
  d->printStepPointMCCollection("g4run","tracker");
}

//_____________________________________________________________________________
void print_Tracks() {
  printf("\n");
  TAnaDump* d = TAnaDump::Instance();
  printf("---- TrkPatRec:\n"); 
  d->printKalRepCollection("TrkPatRec","","",0); 
  printf("---- CalPatRec:\n"); 
  d->printKalRepCollection("CalPatRec","","",0);
}

// //_____________________________________________________________________________
// void print_clusters() {
//   TStnNode*  n;
//   TStnEvent* ev;
//   ev = g.x->GetEvent();
//   n  = (TStnNode*) ev->GetListOfInputNodes()->FindObject("ClusterBlock");
//   printf("\n");
//   n->GetDataBlock()->Print("/ces-only");
// }

// //_____________________________________________________________________________
// void print_electrons() {
//   TStnNode*  n;
//   TStnEvent* ev;
//   ev = g.x->GetEvent();
//   n  = (TStnNode*) ev->GetListOfInputNodes()->FindObject("ElectronBlock");
//   printf("\n");
//   n->GetDataBlock()->Print("");
// }

// //_____________________________________________________________________________
// void print_jets() {
//   printf("\n");
//   m_dbg->fJetBlock->Print();
// }

// //_____________________________________________________________________________
// void print_muons() {
//   printf("\n");
//   m_dbg->fMuonBlock->Print();
// }

// //_____________________________________________________________________________
// void print_obsp() {
//   printf("\n");
//   m_dbg->fObspBlock->Print();
// }

// //_____________________________________________________________________________
// void print_genp() {
//   printf("\n");
//   m_dbg->fGenpBlock->Print();
// }

// //_____________________________________________________________________________
// void print_conversions() {
//   printf("\n");
//   m_dbg->fConversionBlock->Print();
// }

// //_____________________________________________________________________________
// void print_tracks() {
//   printf("\n");

//   Int_t xchan, zchan;
//   Double_t *chan_to_x_map = TCesChamber::fgXStrip;
//   Double_t *chan_to_z_map = TCesChamber::fgXWire;
//   Double_t xcoord, zcoord;
  

//   m_dbg->fTrackBlock->Print();
  
// //    for (int trk=0; trk<block->fNTracks; trk++) {

// //      TStnTrack* t = block->Track(trk);

// //      if(t->fXCes[0] < 999.00 && t->fXCes[1] < 999.00) { 

// //        xcoord = t->fXCes[0];
// //        zcoord = t->fXCes[1];
  
// //        Double_t dxmin, dzmin, dx, dz;

// //        dxmin = 1000;
// //        xchan = -1;

// //        for (int i=0; i<=128; i++) {
// //  	dx = TMath::Abs(chan_to_x_map[i]-xcoord);
// //  	if (dx < dxmin) {
// //  	  dxmin   = dx; 
// //  	  xchan = i;
// //  	}
// //  	else break;
// //        }
      
// //        dzmin = 1000;
// //        zchan = -1;

// //        for (int i=0; i<=128; i++){
// //  	dz = TMath::Abs(chan_to_z_map[i]-zcoord);
// //  	if (dz < dzmin) {
// //  	  dzmin   = dz;
// //  	  zchan = i;
// //  	}
// //  	else break;
// //        }
      
// //        printf("SIDE: %1i  WEDGE:%3i  Z: %6.2f  X: %6.2f ==> XCHAN: %3i ZCHAN: %3i\n",
// //  	     t->fSide,
// //  	     t->fWedge,
// //  	     t->fXCes[0],
// //  	     t->fXCes[1],
// //  	     xchan,
// //  	     zchan);
// //      }    
// //    }
// }


// //_____________________________________________________________________________
// void x(int Run, int Event) {
//   g.x->ProcessEvent(Run,Event); 

//   TObjArray* lc = m_cnv->fListOfConversions;

//   cnv = (TStnConversion*) lc->At(0); 

//   trkb = m_dbg->fTrackBlock;
//   pib  = m_dbg->fPi0Block;
//   taub = m_dbg->fTauBlock;
//   jb   = m_dbg->fJetBlock;
//   eb   = m_dbg->fElectronBlock;
//   mb   = m_dbg->fMuonBlock;

//   int nc = lc->GetEntries();
//   printf( "nconv: %i\n" , nc);
//   if (nc > 0) printf(" tracks: %i %i\n", cnv->TrackNumber(0), cnv->TrackNumber(1));
//   if (m_dbg->GetHeaderBlock()->McFlag()) {
//     m_dbg->fObspBlock->Print(15) ; 
//   }

//   TStntuple::PrintTauBlock(g.x->GetEvent(),3);
// }


// //_____________________________________________________________________________
// void tau_mass_04() {

//   TLorentzVector  *tau_mom, *tmom;

//   TLorentzVector  sum;

//   int nax, nst, n04;

//   trkb = m_dbg->fTrackBlock;
//   pib  = m_dbg->fPi0Block;
//   taub = m_dbg->fTauBlock;

//   int ntrk = trkb->NTracks();

//   for (int i=0; i<taub->NTaus(); i++) {
//     TStnTau* tau = taub->Tau(i);
//     TLorentzVector* tau_mom = tau->TrackMomentum();

//     float vz = tau->Zv();
// 				// loop over all the tracks 
//     sum.SetXYZT(0,0,0,0);
//     n04 = 0;
//     for (int it=0; it<ntrk; it++) {
//       TStnTrack* trk = trkb->Track(it);
//       tmom = trk->Momentum();

//       nax = trk->NCotAxSeg(6) ;
//       nst = trk->NCotStSeg(6) ;
//       if ((nax >= 2) && (nst >= 2) && (tmom->Pt() > 0.4)) {
// 	if (TMath::Abs(trk->Z0()-vz) < 5) {

// 	  if (tau_mom->Angle(tmom->Vect()) < TMath::Pi()/10) {
// 	    sum += *tmom;
// 	    n04++;
// 	  }
// 	}
//       }
//     }

//     printf(" itau, n04, mass: %2i %2i %10.3f\n",i,n04,sum.M());
//   }
// }


// //_____________________________________________________________________________
// void print_tracks_30deg() {

//   TLorentzVector  *tau_mom, *tmom;

//   TLorentzVector  sum;

//   trkb = m_dbg->fTrackBlock;
//   pib  = m_dbg->fPi0Block;
//   taub = m_dbg->fTauBlock;

//   int nax, nst;
//   int ntrk = trkb->NTracks();

//   printf("\n");

//   for (int i=0; i<taub->NTaus(); i++) {
//     TStnTau* tau = taub->Tau(i);
//     TLorentzVector* tau_mom = tau->TrackMomentum();
//     int banner_printed = 0;
// 				// loop over all the tracks 

//     for (int it=0; it<ntrk; it++) {
//       TStnTrack* trk = trkb->Track(it);
//       tmom = trk->Momentum();

//       nax = trk->NCotAxSeg(6) ;
//       nst = trk->NCotStSeg(6) ;
//       if ((nax >= 2) && (nst >= 2)) {
// 	if (tau_mom->Angle(tmom->Vect()) < TMath::Pi()/6.) {
// 	  if (! banner_printed) {
// 	    trk->Print("banner");
// 	    banner_printed=1;
// 	  }
// 	  trk->Print("data");
// 	}
//       }
//     }
//   }
// }


/////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________________
void help() {
//
// This macro generates a Controlbar menu: 
// To see the output, click begin_html <a href="gif/demos.gif" >here</a> end_html
// To execute an item, click with the left mouse button.
// To see the HELP of a button, click on the right mouse button.

   gStyle->SetScreenFactor(1); //if you have a large screen, select 1,2 or 1.4

   bar = new TControlBar("vertical", "Help",10,10);

   //   bar->AddButton("next event", "next_event();","next event");
   bar->AddButton("print StrawHits"   , "print_StrawHitCollection();"   ,"print straw hits"  );
   bar->AddButton("print StepPointMCs", "print_StepPointMCCollection();","print StepPointMCs");
   bar->AddButton("print tracks"      , "print_Tracks();"               ,"print Tracks"      );
//    bar->AddButton("print jets"     , "print_jets();"     ,"print jets");
//    bar->AddButton("print taus",  "print_taus()", "print taus");
//    bar->AddButton("print clusters", "print_clusters()", "print clusters");
//    bar->AddButton("print tracks", "print_tracks()", "print tracks");
//    bar->AddButton("tau mass(04)", "tau_mass_04()", "tau mass");
//    bar->AddButton("print tracks 30 deg", "print_tracks_30deg()", "tracks 30 deg");
//    bar->AddButton("print OBSP", "print_obsp()", "print OBSP block");
//    bar->AddButton("print GENP", "print_genp()", "print GENP block");
//    bar->AddButton("print conversions", "print_conversions()", "print conversion block");
   bar->Show();
   gROOT->SaveContext();
}

