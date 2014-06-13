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

#include "Stntuple/obj/TDisk.hh"
#include "Stntuple/obj/TStnCrystal.hh"

// #include "CalorimeterGeom/inc/VaneCalorimeter.hh"
// #include "CalorimeterGeom/inc/Crystal.hh"
// #include "CalorimeterGeom/inc/Disk.hh"
// #include "CalorimeterGeom/inc/DiskCalorimeter.hh"
// #include "CalorimeterGeom/inc/Calorimeter.hh"

ClassImp(TDisk)

//_____________________________________________________________________________
TDisk::TDisk(): TObject() {
}

//_____________________________________________________________________________
TDisk::TDisk(int SectionID, double RMin, double RMax, double Z0, 
	     double HexSize, double DeadSpace, double MinFraction) : TObject() {

  TStnCrystal     *crystal;
  int             n1;
  TVector2        pos;
  
  double          size;

  fSectionID    = SectionID;
//------------------------------------------------------------------------------
// define number of rings and the total number of crystals with a safety margin, 
// count only crystals fully inside
//-----------------------------------------------------------------------------
  fRMin         = RMin;
  fRMax         = RMax;
  fZ0           = Z0;
  fHexSize      = HexSize;		// full size
  fDeadSpace    = DeadSpace;
  fMinFraction  = MinFraction;
					// account for spacing between the crystals
  size          = fHexSize+2*fDeadSpace;

  fNRings       = (int(fRMax/size-0.2)+1)*1.2+1;

  n1            = 1 + 3*fNRings*(fNRings-1);

  fFirst            = new int[fNRings];
  fNCrystalsPerRing = new int[fNRings];
  fNInside          = new int[fNRings];

  fFirst           [0] = 0;
  fNCrystalsPerRing[0] = 1;
  fNInside         [0] = 0;

  for (int i=1; i<fNRings; i++) {
    fFirst           [i] = 1 + 3*i*(i-1);
    fNCrystalsPerRing[i] = 6*i;
    fNInside         [i] = 0;
  }

  int        loc, inside; 
  double     fraction;
  THexIndex  hex_index;

  fListOfCrystals = new TObjArray(n1);
  fNCrystals   = 0;

  for (int ir=0; ir<fNRings; ir++) {

    for (int ic=0; ic<fNCrystalsPerRing[ir]; ic++) {
      loc       = fFirst[ir]+ic;
      GetHexIndex(loc,&hex_index);
      //      ring      = GetRing  (&hex_index);
      inside    = IsInside(&hex_index,&fraction);

      if (inside) {
	GetPosition(&hex_index,&pos);
	crystal = new TStnCrystal(&hex_index,pos.X(),pos.Y(),Z0,fHexSize);
	crystal->SetDisk(this);
	fListOfCrystals->Add(crystal);
	fNInside[ir] += 1;
      }
    }

    fNCrystals += fNInside[ir];
  }

  //  fFirst          = SectionID*fNCrystals;
}

//-----------------------------------------------------------------------------
TDisk::~TDisk() {
}

//-----------------------------------------------------------------------------
void TDisk::Clear(Option_t* Opt) {

  TStnCrystal  *cr;
  
  for (int i=0; i<fNCrystals; i++) {
    cr     = (TStnCrystal*) fListOfCrystals->UncheckedAt(i);
    cr->Clear();
  }
}

//-----------------------------------------------------------------------------
void TDisk::Print(Option_t* Opt) const {
  printf("------------------------------------------------------------------------\n");
  printf(" ID   RMin    RMax   HexSize   MaxFraction NCrystals   NRings  Offset  \n");
  printf("------------------------------------------------------------------------\n");
  printf("%3i %8.3f %8.3f %8.3f  %5.1f %5i  %5i  %5i\n",
	 fSectionID,fRMin,fRMax,fHexSize,fMinFraction,fNCrystals,
	 fNRings,fChannelOffset);
}

//-----------------------------------------------------------------------------
void TDisk::GetHexIndex(int I, THexIndex* Index) {

  int      ring, loc, l(0), k(0), n1, n2;

  if (I > 0) {
				// find the ring number
    for (int i=1; i>0; i++) {
      n1 = 3*(i-1)*i+1;
      n2 = 3*i*(i+1)+1;
      if (I < n2) {
	ring = i;
	loc  = I-n1;
	break;
      }
    }

    int seg1 = loc / ring; 
    int pos  = loc % ring;

    int seg2 = seg1+1;
    if (seg2 == 6) seg2 = 0;

    l = THexIndex::fgPos[seg1].fL*ring + (THexIndex::fgPos[seg2].fL-THexIndex::fgPos[seg1].fL)*pos;
    k = THexIndex::fgPos[seg1].fK*ring + (THexIndex::fgPos[seg2].fK-THexIndex::fgPos[seg1].fK)*pos;
  }
  Index->Set(l,k);
}

//-----------------------------------------------------------------------------
int TDisk::GetRing(THexIndex* Index) {

  int  ring;

  if (Index->fL*Index->fK > 0) {
    ring = std::abs(Index->fL+Index->fK);
  }
  else if ( std::abs(Index->fL) > std::abs(Index->fK) ) {
    ring = std::abs(Index->fL);
  }
  else {
    ring = std::abs(Index->fK);
  }

  return ring;
}

//-----------------------------------------------------------------------------
int TDisk::GetRing(int I) {
  THexIndex  hex_index;
  
  GetHexIndex(I,&hex_index);

  return GetRing(&hex_index);
}

//-----------------------------------------------------------------------------
void TDisk::GetPosition(THexIndex* Index, TVector2* Pos) {
  double x, y, step;

  step = fHexSize+2*fDeadSpace;
  x    = step*(Index->fL+Index->fK)*sqrt(3.)/2.;
  y    = step*(Index->fL-Index->fK)/2.;
  Pos->Set(x,y);

}

//-----------------------------------------------------------------------------
void TDisk::GetPosition(int I, TVector2* Pos) {
  THexIndex  hex_index;
  
  GetHexIndex(I,&hex_index);
  GetPosition(&hex_index,Pos);
}

//-----------------------------------------------------------------------------
int TDisk::IsInside(THexIndex* Index, double* Fraction) {

  int       inside, nvin(0), nbelow(0), nabove(0);
  double    x0, y0, r, x, y, phi, s, s0, s1, r0, dr, adr;
  double    fr, step;
					// crystal center position
  step = fHexSize+2*fDeadSpace;
  x0   = step*(Index->fL+Index->fK)*sqrt(3.)/2.;
  y0   = step*(Index->fL-Index->fK)/2.;

					// loop over 6 vertices and check if 
					// all of them are inside
  double rho = fHexSize/sqrt(3.);

  for (int i=0; i<6; i++) {
    phi = i*TMath::Pi()/3;
    x   = x0+rho*TMath::Cos(phi);
    y   = y0+rho*TMath::Sin(phi);

    r   = sqrt(x*x+y*y);

    if (r < fRMin) {
      nbelow += 1;
    }
    else if (r > fRMax) {
      nabove += 1;
    }
    else {
      nvin += 1;
    }
  }
//-----------------------------------------------------------------------------
// 'nvin' - number of vertices inside the disk annulus
//-----------------------------------------------------------------------------
  if (nvin == 0) {
    fr = 0;
  }
  else if (nvin == 6) {
    fr = 1.;
  }
  else { 
//-----------------------------------------------------------------------------
// at this point can calculate fraction of the crystal area inside the ring
//-----------------------------------------------------------------------------
    r0  = sqrt(x0*x0+y0*y0);
    dr  = r0-fRMin;

    if (nabove > 0) {
//-----------------------------------------------------------------------------
// crystal crosses the outer ring
//-----------------------------------------------------------------------------
      dr  = r0 -fRMax;
      adr = fabs(dr);
      if (adr > fHexSize/2) adr = fHexSize/2.;

      s   = fHexSize*fHexSize*sqrt(3)/2;
      s1  = (2*fHexSize-adr)*adr/sqrt(3.);
      s0  = (3*fHexSize-2*adr)*(fHexSize-2*adr)/4/sqrt(3);

      if (dr <= 0) {
	fr = 1-s0/s;
      }
      else {
        fr = s0/s;
      }
    }
    else {
//-----------------------------------------------------------------------------
// crystal crosses the inner ring
//-----------------------------------------------------------------------------
      adr = fabs(dr);
      if (adr > fHexSize/2) adr = fHexSize/2.;

      s   = fHexSize*fHexSize*sqrt(3)/2;
      s1  = (2*fHexSize-adr)*adr/sqrt(3.);
      s0  = (3*fHexSize-2*adr)*(fHexSize-2*adr)/4/sqrt(3);
      if (dr > 0) {
	fr = 1-s0/s;
      }
      else {
        fr = s0/s;
      }
    }
  }

  if (fr >= fMinFraction) inside = 1;
  else                    inside = 0;

  *Fraction = fr;

  return inside;
}

//-----------------------------------------------------------------------------
void TDisk::Paint(Option_t* Opt) {

  TVector2 pos;

  THexIndex*  hex_index;

  TPolyLine p;
  TEllipse  e;

  double x[7], y[7];

  //  TCanvas* c = new TCanvas("c","c",800,800);

  // fListOfPolyLines->Delete();

  TStnCrystal* cr;

  for (int i=0; i<fNCrystals; i++) {
    cr = (TStnCrystal*) fListOfCrystals->At(i);
    hex_index = cr->HexIndex();

    pos.Set(cr->Center()->X(),cr->Center()->Y());
    double   x0, y0;

    x0 = cr->Center()->X();
    y0 = cr->Center()->Y();

    //    GetPosition(i,&pos);

    //    int ring = GetRing(i);

    //    printf(" i, r, l,k, x,y : %5i %5i %5i %5i %10.4lf %10.4lf\n",i,ring, hex_index.fL,hex_index.fK,pos.X(),pos.Y());

    double rho = fHexSize/sqrt(3.);

    for (int iv=0; iv<7; iv++) {
      int l = iv % 6;
      x[iv] = x0 + rho*cos(TMath::Pi()*l/3);
      y[iv] = y0 + rho*sin(TMath::Pi()*l/3);
    }

    p.SetLineColor(1);
    p.SetLineWidth(1);

    double fr;

    if (IsInside(hex_index,&fr)) {

      if (fr == 1) {
	p.SetFillColor(2);                     // fully inside
	p.SetFillStyle(3001);
      }
      else {
	p.SetFillColor(kBlue);             // partially inside
	p.SetFillStyle(1001);
      }

      //      p.SetFillStyle(3001);
      p.DrawPolyLine(7,x,y,"F");
    }
    else {
				// fully outside
      p.SetFillStyle(0);
      p.SetFillColor(0);
    }

    p.DrawPolyLine(7,x,y,"");
  } 

  e.SetFillStyle(0);
  e.SetFillColor(0);
  e.DrawEllipse(0,0,fRMin,fRMin,0,360,0,"");
  e.DrawEllipse(0,0,fRMax,fRMax,0,360,0,"");
  

//   c->Update();
//   c->Draw();
}
