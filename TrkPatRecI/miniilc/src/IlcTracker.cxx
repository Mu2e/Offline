/**************************************************************************
 * Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 *                                                                        *
// Author: The ILC Off-line Project. 
 // Part of the code has been developed by Alice Off-line Project. 
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: IlcTracker.cxx,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//-------------------------------------------------------------------------
//               Implementation of the IlcTracker class
//  that is the base for IlcTPCtracker, IlcVXDtrackerV2 and IlcTRDtracker    
//        Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------
#include <TClass.h>
#include <TMath.h>


//#include "IlcMagF.h"
#include "IlcTracker.h"
#include "IlcCluster.h"
#include "IlcKalmanTrack.h"

#include <iostream>
#include <map>


Bool_t IlcTracker::fgUniformField=kTRUE;
Double_t IlcTracker::fgBz=kAlmost0Field;
//const IlcMagF *IlcTracker::fgkFieldMap=0;

//ClassImp(IlcTracker)

IlcTracker::IlcTracker():
  TObject(),
  fX(0),
  fY(0),
  fZ(0),
  fSigmaX(0.005),
  fSigmaY(0.005),
  fSigmaZ(0.010)
{
  //--------------------------------------------------------------------
  // The default constructor.
  //--------------------------------------------------------------------
  //  if (!fgkFieldMap) IlcWarning("Field map is not set. Call IlcTracker::SetFieldMap before creating a tracker!");
}

//__________________________________________________________________________
IlcTracker::IlcTracker(const IlcTracker &atr):
  TObject(atr),
  fX(atr.fX),
  fY(atr.fY),
  fZ(atr.fZ),
  fSigmaX(atr.fSigmaX),
  fSigmaY(atr.fSigmaY),
  fSigmaZ(atr.fSigmaZ)
{
  //--------------------------------------------------------------------
  // The default constructor.
  //--------------------------------------------------------------------
  //  if (!fgkFieldMap) IlcWarning("Field map is not set. Call IlcTracker::SetFieldMap before creating a tracker!");
}
/*
//__________________________________________________________________________
void IlcTracker::SetFieldMap(const IlcMagF* map, Bool_t uni) {
  //--------------------------------------------------------------------
  //This passes the field map to the reconstruction.
  //--------------------------------------------------------------------
  if (map==0) IlcFatalClass("Can't access the field map !");

  if (fgkFieldMap) {
     IlcWarningClass("The magnetic field map has been already set !");
     //     return;
  }

  fgUniformField=uni;
  fgkFieldMap=map;

  //Float_t r[3]={0.,0.,0.},b[3]; map->Field(r,b);
  //Double_t bz=-b[2];
 
  Double_t bz=-map->SolenoidField();
  fgBz=TMath::Sign(1e-13,bz) + bz;
  float x[4]={0,0,0,0},b[4];
  map->Field(x,b);
  IlcDebugGeneral("IlcTracker",1,Form("Field in midpoint B=(%f,%f,%f), uniform=%f",b[0],b[1],b[2],bz));
  if(fabs(b[2]+bz)>1e-6){
    IlcWarningGeneral("IlcTracker",
		      Form("Differnt field in midpoint map=%f, uniform=%f; I will set from midpoint ",b[2],-bz));
    bz=-b[2];
  }
  fgBz=TMath::Sign(kAlmost0Field,bz) + bz;

}
*/
//__________________________________________________________________________
void IlcTracker::CookLabel(IlcKalmanTrack *t, Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  std::map<int,int> labels;
  int numberofpoints=t->GetNumberOfClusters();
  for(int i=0;i<numberofpoints;i++) {
    Int_t cindex = t->GetClusterIndex(i);
    IlcCluster *cl = GetCluster(cindex);
    labels[cl->GetLabel(0)]++;
  }
  int maxlabel=-1,ntimes=0;
  for(std::map<int,int>::iterator it=labels.begin();it!=labels.end();it++){
    if((*it).second>ntimes){
      maxlabel=(*it).first;
      ntimes=(*it).second;
    }
  }
  Int_t lflag=0;
  
  for(Int_t i=0;i<numberofpoints;i++){
    Int_t cindex = t->GetClusterIndex(i);
    IlcCluster *cl = GetCluster(cindex);
    
    if(cl->GetLabel(0)==maxlabel || 
       cl->GetLabel(1)==maxlabel || 
       cl->GetLabel(2)==maxlabel) lflag++;
  }
  
  double fakeratio=1.-double(lflag)/numberofpoints;
  if(maxlabel<0) maxlabel=1000000-maxlabel;

  t->SetFakeRatio(fakeratio);
  t->SetLabel(t->GetFakeRatio()>wrong?-maxlabel:maxlabel);
  
}

//____________________________________________________________________________
void IlcTracker::UseClusters(const IlcKalmanTrack *t, Int_t from) const {
  //------------------------------------------------------------------
  //This function marks clusters associated with the track.
  //------------------------------------------------------------------
  Int_t noc=t->GetNumberOfClusters();
  for (Int_t i=from; i<noc; i++) {
     Int_t index=t->GetClusterIndex(i);
     IlcCluster *c=GetCluster(index); 
     c->Use();   
  }
}

Double_t IlcTracker::GetBz(Float_t *r) {
  //------------------------------------------------------------------
  // Returns Bz (kG) at the point "r" .
  //------------------------------------------------------------------
  return GetBz();
  //    Float_t b[3]; fgkFieldMap->Field(r,b);
  //  Double_t bz=-Double_t(b[2]);
  //  return  (TMath::Sign(kAlmost0Field,bz) + bz);
}
