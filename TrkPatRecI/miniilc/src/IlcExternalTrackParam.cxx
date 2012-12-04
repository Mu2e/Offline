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

/* $Id: IlcExternalTrackParam.cxx,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Implementation of the external track parameterisation class.              //
//                                                                           //
// This parameterisation is used to exchange tracks between the detectors.   //
// A set of functions returning the position and the momentum of tracks      //
// in the global coordinate system as well as the track impact parameters    //
// are implemented.
// Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                            //
///////////////////////////////////////////////////////////////////////////////
#include "IlcExternalTrackParam.h"
//#include "IlcESDVertex.h"
#include "IlcLog.h"
#include "IlcKalmanTrack.h"
#include <TMatrixDSym.h>
#include <TMath.h>
#include <math.h>
#include <TClass.h>
#include <iostream>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
using namespace std; 

//ClassImp(IlcExternalTrackParam)

//_____________________________________________________________________________
IlcExternalTrackParam::IlcExternalTrackParam() :
  TObject(),
  kEndCap(false),
  fX(0),
  fAlpha(0)
{
  //
  // default constructor
  fDir=1.;
  //
  for (Int_t i = 0; i < 5; i++) fP[i] = 0;
  for (Int_t i = 0; i < 15; i++) fC[i] = 0;
}

//_____________________________________________________________________________
IlcExternalTrackParam::IlcExternalTrackParam(const IlcExternalTrackParam &track):
  TObject(track),
  kEndCap(track.kEndCap),
  fX(track.fX),
  fAlpha(track.fAlpha)
{
  //
  // copy constructor
  //
  for (Int_t i = 0; i < 5; i++) fP[i] = track.fP[i];
  fDir=track.fDir;
  if(fabs(fDir)!=1.) fDir=TMath::Sign(1.,fDir);
  for (Int_t i = 0; i < 15; i++) fC[i] = track.fC[i];
}

//_____________________________________________________________________________
IlcExternalTrackParam::IlcExternalTrackParam(Double_t x, Double_t alpha, 
					     const Double_t param[5], 
					     const Double_t covar[15],
					     double fDirection) :
  TObject(),
  kEndCap(false),
  fX(x),
  fAlpha(alpha)
{
  //
  // create external track parameters from given arguments
  fDir=fDirection;
  //
  for (Int_t i = 0; i < 5; i++)  fP[i] = param[i];
  for (Int_t i = 0; i < 15; i++) fC[i] = covar[i];
}

//_____________________________________________________________________________
void IlcExternalTrackParam::Set(Double_t x, Double_t alpha,
				const Double_t p[5], const Double_t cov[15],
				double fDirection) {
  //
  //  Sets the parameters
  //
  fX=x;
  fAlpha=alpha;
  for (Int_t i = 0; i < 5; i++)  fP[i] = p[i];
  if(fabs(fP[2])>kAlmost1) fP[2]=TMath::Sign(kAlmost1,fP[2]);
  fDir=fDirection;
  for (Int_t i = 0; i < 15; i++) fC[i] = cov[i];
}

void IlcExternalTrackParam::Set(const IlcExternalTrackParam *t){
  IlcExternalTrackParam::operator=(*t);
  if(IlcDebugLevelClass()>5){
    IlcDebug(6,"Set IlxExternal Params to valuse:")
    IlcExternalTrackParam::Print("param");
  }
}

//_____________________________________________________________________________
void IlcExternalTrackParam::Reset() {
  //
  // Resets all the parameters to 0 
  //
  fX=fAlpha=0.;
  for (Int_t i = 0; i < 5; i++) fP[i] = 0;
  fDir=1.;
  for (Int_t i = 0; i < 15; i++) fC[i] = 0;
}

Double_t IlcExternalTrackParam::GetP() const {
  //---------------------------------------------------------------------
  // This function returns the track momentum
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  if (TMath::Abs(fP[4])<=kAlmost0) return kVeryBig;
  return TMath::Sqrt(1.+ fP[3]*fP[3])/TMath::Abs(fP[4]);
}

Double_t IlcExternalTrackParam::Get1P() const {
  //---------------------------------------------------------------------
  // This function returns the 1/(track momentum)
  //---------------------------------------------------------------------
  return TMath::Abs(fP[4])/TMath::Sqrt(1.+ fP[3]*fP[3]);
}

//_______________________________________________________________________
Double_t IlcExternalTrackParam::GetD(Double_t x,Double_t y,Double_t b) const {
  //------------------------------------------------------------------
  // This function calculates the transverse impact parameter
  // with respect to a point with global coordinates (x,y)
  // in the magnetic field "b" (kG)
  //------------------------------------------------------------------
  if (TMath::Abs(b) < kAlmost0Field) return GetLinearD(x,y);
  Double_t rp4=GetC(b);

  Double_t xt=fX, yt=fP[0];

  Double_t sn=TMath::Sin(fAlpha), cs=TMath::Cos(fAlpha);
  Double_t a = x*cs + y*sn;
  y = -x*sn + y*cs; x=a;
  xt-=x; yt-=y;

  sn=rp4*xt - fP[2]; cs=rp4*yt + fDir*TMath::Sqrt(1.- fP[2]*fP[2]);
  a=2*(xt*fP[2] - yt*fDir*TMath::Sqrt(1.- fP[2]*fP[2]))-rp4*(xt*xt + yt*yt);
  return  -a/(1 + TMath::Sqrt(sn*sn + cs*cs));
}

//_______________________________________________________________________
void IlcExternalTrackParam::
GetDZ(Double_t x, Double_t y, Double_t z, Double_t b, Float_t dz[2]) const {
  //------------------------------------------------------------------
  // This function calculates the transverse and longitudinal impact parameters
  // with respect to a point with global coordinates (x,y)
  // in the magnetic field "b" (kG)
  //------------------------------------------------------------------
  Double_t f1 = fP[2], r1 = fDir*TMath::Sqrt(1. - f1*f1);
  Double_t xt=fX, yt=fP[0];
  Double_t sn=TMath::Sin(fAlpha), cs=TMath::Cos(fAlpha);
  Double_t a = x*cs + y*sn;
  y = -x*sn + y*cs; x=a;
  xt-=x; yt-=y;

  Double_t rp4=GetC(b);
  if ((TMath::Abs(b) < kAlmost0Field) || (TMath::Abs(rp4) < kAlmost0)) {
     dz[0] = -(xt*f1 - yt*r1);
     dz[1] = fP[1] + (dz[0]*f1 - xt)/r1*fP[3] - z;
     return;
  }

  sn=rp4*xt - f1; cs=rp4*yt + r1;
  a=2*(xt*f1 - yt*r1)-rp4*(xt*xt + yt*yt);
  Double_t rr=TMath::Sqrt(sn*sn + cs*cs);
  dz[0] = -a/(1 + rr);
  Double_t f2 = -sn/rr, r2 = fDir*TMath::Sqrt(1. - f2*f2);
  dz[1] = fP[1] + fP[3]/rp4*TMath::ASin(f2*r1 - f1*r2) - z;
}

//_______________________________________________________________________
Double_t IlcExternalTrackParam::GetLinearD(Double_t xv,Double_t yv) const {
  //------------------------------------------------------------------
  // This function calculates the transverse impact parameter
  // with respect to a point with global coordinates (xv,yv)
  // neglecting the track curvature.
  //------------------------------------------------------------------
  Double_t sn=TMath::Sin(fAlpha), cs=TMath::Cos(fAlpha);
  Double_t x= xv*cs + yv*sn;
  Double_t y=-xv*sn + yv*cs;

  Double_t d = (fX-x)*fP[2] - (fP[0]-y)*TMath::Sqrt(1.- fP[2]*fP[2])*fDir;

  return -d;
}

Bool_t IlcExternalTrackParam::CorrectForMaterial
(Double_t d,  Double_t x0, Double_t mass, Double_t (*Bethe)(Double_t)) {
  //------------------------------------------------------------------
  // This function corrects the track parameters for the crossed material
  // "d"    - the thickness (fraction of the radiation length)
  // "x0"   - the radiation length (g/cm^2) 
  // "mass" - the mass of this particle (GeV/c^2)
  //------------------------------------------------------------------
  Double_t &fP2=fP[2];
  Double_t &fP3=fP[3];
  Double_t &fP4=fP[4];

  Double_t &fC22=fC[5];
  Double_t &fC33=fC[9];
  Double_t &fC43=fC[13];
  Double_t &fC44=fC[14];

  Double_t p=GetP();
  Double_t p2=p*p;
  Double_t beta2=1.-mass*mass/(p2 + mass*mass);

  if(fabs(GetP())<kSmallPt) return kFALSE; //cut on 1MeV        

  //Multiple scattering******************
  if (d!=0) {
     Double_t theta2=14.1*14.1/(beta2*p2*1e6)*TMath::Abs(d);
     //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*TMath::Abs(d)*9.36*2.33;
     fC22 += theta2*(1.- fP2*fP2)*(1. + fP3*fP3);
     fC33 += theta2*(1. + fP3*fP3)*(1. + fP3*fP3);
     fC43 += theta2*fP3*fP4*(1. + fP3*fP3);
     fC44 += theta2*fP3*fP4*fP3*fP4;
  }

  //Energy losses************************
  if (x0!=0. && beta2<1.) {
     d*=x0;
     Double_t dE=Bethe(beta2)*d;
     Double_t e=TMath::Sqrt(p2 + mass*mass);
     if(dE*e>=p2) return kFALSE;
     fP4*=(1.- e/p2*dE);
     if(fabs(GetP())<kSmallPt) return kFALSE; //cut on 1MeV

     // Approximate energy loss fluctuation (M.Ivanov)
     const Double_t cnst=0.07; // To be tuned.  
     Double_t sigmadE=cnst*TMath::Sqrt(TMath::Abs(dE)); 
     fC44+=((sigmadE*e/p2*fP4)*(sigmadE*e/p2*fP4)); 
 
  }

  return kTRUE;
}

Double_t ApproximateBetheBloch(Double_t beta2) {
  //------------------------------------------------------------------
  // This is an approximation of the Bethe-Bloch formula with 
  // the density effect taken into account at beta*gamma > 3.5
  // (the approximation is reasonable only for solid materials) 
  //------------------------------------------------------------------
  if(beta2>0.9999999999999) beta2=0.9999999999999;
  if (beta2/(1-beta2)>3.5*3.5)
     return 0.153e-3/beta2*(log(3.5*5940)+0.5*log(beta2/(1-beta2)) - beta2);

  return 0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2);
}

bool  IlcExternalTrackParam::ChangeDirection(){
  // change of track direction 
  fDir*=-1;
  fP[2]*=-1;
  fP[3]*=-1;
  fP[4]*=-1;
  
  for(int i=2;i<5;i++){
    for(int j=0;j<5;j++){
      if(j!=i) Cov(i,j)*=-1;
    }
  }
  return true;
};


Bool_t IlcExternalTrackParam::Rotate(Double_t alpha,Double_t bz) {
  //------------------------------------------------------------------
  // Transform this track to the local coord. system rotated
  // by angle "alpha" (rad) with respect to the global coord. system. 
  //------------------------------------------------------------------
  if (TMath::Abs(fP[2]) >= kAlmost1) {
    //     IlcError(Form("Precondition is not satisfied: |sin(phi)|>1 ! %f",fP[2])); 
     return kFALSE;
  }
  if      (alpha < -TMath::Pi()) alpha += 2*TMath::Pi();
  else if (alpha >= TMath::Pi()) alpha -= 2*TMath::Pi();

  Double_t &fP0=fP[0];
  Double_t &fP2=fP[2];
  Double_t &fP3=fP[3];
  Double_t &fC00=fC[0];
  Double_t &fC10=fC[1];
  Double_t &fC20=fC[3];
  Double_t &fC30=fC[6];
  Double_t &fC40=fC[10];
  Double_t &fC11=fC[2];
  Double_t &fC21=fC[4];
  Double_t &fC31=fC[7];
  Double_t &fC41=fC[11];
  Double_t &fC22=fC[5];
  Double_t &fC32=fC[8];
  Double_t &fC42=fC[12];

  Double_t x=fX;
  Double_t ca,sa;sincos(alpha-fAlpha,&sa,&ca);
  Double_t sf=fP2, cf=fDir*TMath::Sqrt(TMath::Abs(1.- fP2*fP2));

  Double_t tmp=sf*ca - cf*sa;
  if (TMath::Abs(tmp) >= kAlmost1) return kFALSE;
  if(fDir*(cf*ca+sf*sa)<0){
    fDir*=-1;
    //    IlcDebug(5,Form("Warning change of direction before rotation cos(direction)= %f after =%f ",cf,(cf*ca+sf*sa)));
  }

  fAlpha = alpha;
  fX =  x*ca + fP0*sa;
  fP0= -x*sa + fP0*ca;
  fP2=  tmp;

  if (TMath::Abs(cf)<kAlmost0) {
    //    IlcError(Form("Too small cosine value %f",cf)); 
    cf = kAlmost0;
  } 

  Double_t newcf=fDir*TMath::Sqrt(TMath::Abs(1.- fP2*fP2));
  if (TMath::Abs(newcf)<kAlmost0) {
    //    IlcError(Form("Too small new cosine value %f",newcf)); 
    newcf = kAlmost0;
  } 
  
  Double_t crv=GetC(bz);

  // F=dg/dp+df(g)/dx*dX0/dp
  //f = F - 1
  Double_t f00=ca - fP2/newcf*sa - 1, f10 = - fP3/newcf*sa , f20 = -crv*sa,   
           f22=(ca+sf/cf*sa)-1;

  //b = C*ft
  Double_t b00=fC00*f00, b01=fC00*f10, b02=fC00*f20+fC20*f22;
  Double_t b10=fC10*f00, b11=fC10*f10, b12=fC10*f20+fC21*f22;
  Double_t b20=fC20*f00, b21=fC20*f10, b22=fC20*f20+fC22*f22;
  Double_t b30=fC30*f00, b31=fC30*f10, b32=fC30*f20+fC32*f22;
  Double_t b40=fC40*f00, b41=fC40*f10, b42=fC40*f20+fC42*f22;

  //a = f*b = f*C*ft
  Double_t a00=f00*b00, a01 = f00*b01, a02=f00*b02, 
                        a11 = f10*b01, a12=f10*b02,
                                       a22=f20*b02+f22*b22;

  //F*C*Ft = C + (a + b + bt)
  fC00 += a00 + 2*b00;
  fC10 += a01 + b10 + b01;
  fC20 += a02+b20+b02;
  fC30 += b30;
  fC40 += b40;
  fC11 += a11 + 2*b11;
  fC21 += a12+ b12+b21;
  fC31 += b31;
  fC41 += b41;
  fC22 += a22 + 2*b22;
  fC32 += b32;
  fC42 += b42;
  return kTRUE;

  /*  Double_t rr=(ca+sf/cf*sa);  

  fC00 *= (ca*ca);
  fC10 *= ca;
  fC20 *= ca*rr;
  fC21 *= rr;
  fC22 *= rr*rr;
  fC30 *= ca;
  fC32 *= rr;
  fC40 *= ca;
  fC42 *= rr;

  return kTRUE;*/
}

Bool_t IlcExternalTrackParam::RotateZ(bool toendcap,Double_t bz) {
  //------------------------------------------------------------------
  // Transform this track to the local coord. system rotated
  // by angle "alpha" (rad) with respect to the global coord. system. 
  //------------------------------------------------------------------
  if(toendcap==kEndCap){
    //     IlcError(Form("Already rotated in Z : %i ",int(kEndCap))); 
     return kFALSE;    
  }

  if (TMath::Abs(fP[2]) >= kAlmost1) {
    //     IlcError(Form("Precondition is not satisfied: |sin(phi)|>1 ! %f",fP[2])); 
     return kFALSE;
  }

  Double_t &fP2=fP[2];
  Double_t &fP3=fP[3];
  Double_t &fC00=fC[0];
  Double_t &fC10=fC[1];
  Double_t &fC20=fC[3];
  Double_t &fC30=fC[6];
  Double_t &fC40=fC[10];
  Double_t &fC11=fC[2];
  Double_t &fC21=fC[4];
  Double_t &fC31=fC[7];
  Double_t &fC41=fC[11];
  Double_t &fC22=fC[5];
  Double_t &fC32=fC[8];
  Double_t &fC42=fC[12];

  Double_t sf=fP2, cf=fDir*TMath::Sqrt(TMath::Abs(1.- fP2*fP2));


  if (TMath::Abs(fP3)<kAlmost0) {
    //    IlcError(Form("Too small tangent %f",fP3)); 
    return false;
  } 
  if (TMath::Abs(cf)<kAlmost0) {
    //    IlcError(Form("Too small cosine %f",cf)); 
    cf = kAlmost0;
  } 

  Double_t crv=GetC(bz);

  // F=dg/dp+df(g)/dx*dX0/dp
  //f = F - 1
  Double_t f01,f11,f21;
  if(toendcap){
    f01=-sf/fP3, f11 =-1 - cf/fP3, f21 = -cf/fP3*crv;
    kEndCap=true;
  }else{
    f01=-sf/cf, f11 =-1 - 1./cf*fP3, f21 = -crv;
    kEndCap=false;
  }
  //b = C*ft
  Double_t b00=fC10*f01, b01=fC10*f11, b02=fC10*f21;
  Double_t b10=fC11*f01, b11=fC11*f11, b12=fC11*f21;
  Double_t b20=fC21*f01, b21=fC21*f11, b22=fC21*f21;
  Double_t b30=fC31*f01, b31=fC31*f11, b32=fC31*f21;
  Double_t b40=fC41*f01, b41=fC41*f11, b42=fC41*f21;

  //a = f*b = f*C*ft
  Double_t a00=f01*b10, a01 = f01*b11, a02=f01*b12, 
                        a11 = f11*b11, a12=f11*b12,
                                       a22=f21*b12;

  //F*C*Ft = C + (a + b + bt)
  fC00 += a00 + 2*b00;
  fC10 += a01 + b10 + b01;
  fC20 += a02+b20+b02;
  fC30 += b30;
  fC40 += b40;
  fC11 += a11 + 2*b11;
  fC21 += a12+ b12+b21;
  fC31 += b31;
  fC41 += b41;
  fC22 += a22 + 2*b22;
  fC32 += b32;
  fC42 += b42;
  return kTRUE;

}

double IlcExternalTrackParam::GetNLoopsToZ(double z,double b){
  Double_t crv=GetC(b);
  if (TMath::Abs(b) < kAlmost0Field) crv=0.; 
  double dz=z-fP[1];
  if (TMath::Abs(dz)<=kAlmost0)  return 0.;                                                                
  if (TMath::Abs(fP[3]) <= kAlmost0) return 1e10;                                                           
  Double_t dl=dz/fP[3];                                                                                       
  return crv*dl/TMath::Pi()/2;
}

Bool_t IlcExternalTrackParam::PropagateTo(Double_t xk, Double_t b,bool byR,bool checkSignX) {
  //----------------------------------------------------------------
  // Propagate this track to the plane X=xk (cm) in the field "b" (kG)
  // if byR=true than propagate to xk according to direction along radius
  // (it must be rotatet before to final phi)
  //----------------------------------------------------------------
  double dx,dy,dz;
  Double_t crv=GetC(b);
  if (TMath::Abs(b) < kAlmost0Field) crv=0.;
  if(fabs(crv)<kAlmost0Curv) crv=kAlmost0Curv;
  Double_t cc=b*kB2C;
  double R=1./crv;

  Double_t f1=fP[2];
  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  Double_t r1=fDir*TMath::Sqrt(1.- f1*f1);
  double f02=0,f03=0,f04=0,f12=0,f13=0,f14=0,f22=0,f23=0,f24=0;

  double phidir=TMath::ASin(fP[2]);
  if(fDir<0) phidir=TMath::Pi()-phidir;
  if(kEndCap){
    dz=xk-fP[1];
    if (TMath::Abs(dz)<=kAlmost0)  return kTRUE;
    if (TMath::Abs(fP[3]) <= kAlmost0) return kFALSE;
    Double_t dl=dz/fP[3];
    Double_t sindl,cosdl;sincos(0.5*crv*dl,&sindl,&cosdl);
    phidir+=crv*dl;
    
    fDir=((fabs(Phi_mpi_pi(phidir))<TMath::Pi()/2)?1:-1);

    dx=2.*R*sindl*(cosdl*r1-sindl*f1);
    dy=2.*R*sindl*(cosdl*f1+sindl*r1);

    Double_t f2=f1 + crv*dx;
    if (TMath::Abs(f2) >= kAlmost1) return kFALSE;

    //f = F - 1
    f02=1./r1*dx;
    
    f04=-R*dy+R*dl*(cosdl*(cosdl*f1+sindl*r1)+
		    sindl*(-sindl*f1+cosdl*r1));f04*=cc;
    f03=dl*(cosdl*(cosdl*f1+sindl*r1)+
	    sindl*(-sindl*f1+cosdl*r1));f03*=-1./fP[3];

    f12=-1./r1*dy;
    
    f14=-R*dx+R*dl*(cosdl*(cosdl*r1-sindl*f1)+
			      sindl*(-sindl*r1-cosdl*f1));f14*=cc;
    f13=dl*(cosdl*(cosdl*r1-sindl*f1)+
	    sindl*(-sindl*r1-cosdl*f1));f13*=-1./fP[3];
    
    f22=crv*f12;
    f23=crv*f13;
    f24=crv*f14+dx*cc;
  }else{
    dx=xk-fX;
    if (TMath::Abs(dx)<=kAlmost0&&!byR)  return kTRUE;
    Double_t f2=f1 + crv*dx;
    if (TMath::Abs(f2) >= kAlmost1) return kFALSE;
    Double_t r2=fDir*TMath::Sqrt(1.- f2*f2);
    dy=dx*(f1+f2)/(r1+r2);
    double invsqrt=1./TMath::Sqrt(0.5*(1.+r1*r2-f1*f2))/fDir;
    double sinphi=dx*crv/2*invsqrt;
    if(TMath::Abs(sinphi)>= kAlmost1) return kFALSE;
    double dphi=TMath::ASin(sinphi);
    double invcosphi=1./TMath::Sqrt(1-sinphi*sinphi);
    dz=2.*R*dphi*fP[3];

    //f = F - 1
    f02=    dx/(r1*r2*(r1+r2))*(1.+f1*f2+r1*r2);           
    f04=    dx*dx/(r2*(r1+r2)*(r1+r2))*(1.+f1*f2+r1*r2);  f04*=cc;
    f12=    dx*fP[3]*invcosphi*0.25*(f1*r2+f2*r1)*(r1+r2)/r1/r2*
            invsqrt*invsqrt*invsqrt;
    f14=  -dz*R+ dx*fP[3]*R*invcosphi*0.5*(r1+r2)/r2*invsqrt; f14*=cc;
    f13=    2.*R*dphi;
    f24=    dx;                       f24*=cc;

    if(byR){
      double Rdirection=r1*fX+fP[2]*fP[0];
      double dr=xk-TMath::Hypot(fX,fP[0]);
      if(fabs(fP[0]+dy)>fabs(fP[0]-dy+2*r1*R)){
	dy=-dy+2*r1*R;

	double Ronedirection=(fP[0]+r1*R);
	
	double signdir=(Ronedirection*Rdirection>0)?1.:-1.;
	double signdl=(dr*Rdirection*R>0)?1.:-1.;
	double dl=TMath::ACos(f1*signdir)*2.*R*signdl;
	dz=-dz+dl*fP[3];
	
	f02=-f02-2*f1/r1*R;
	f04=-f04-2*r1*R*R*cc;
	
	f12=-f12-2.*R*signdl*fP[3]/r1*fDir*signdir;
	f13=-f13+dl;
	f14=-f14-R*dl*fP[3]*cc;


	fDir*=-1;
      }
    }
  }



  Double_t &fP0=fP[0], &fP1=fP[1], &fP2=fP[2];//, &fP3=fP[3], &fP4=fP[4];
  Double_t 
  &fC00=fC[0],
  &fC10=fC[1],   &fC11=fC[2],  
  &fC20=fC[3],   &fC21=fC[4],   &fC22=fC[5],
  &fC30=fC[6],   &fC31=fC[7],   &fC32=fC[8],   &fC33=fC[9],  
  &fC40=fC[10],  &fC41=fC[11],  &fC42=fC[12],  &fC43=fC[13], &fC44=fC[14];


  fX+=dx;
  fP0 += dy;
  fP1 += dz;  
  fP2 += dx*crv;

   
  
  //b = C*ft
  Double_t b00=f02*fC20 + f03*fC30 + f04*fC40, b01=f12*fC20 + f14*fC40 + f13*fC30;
  Double_t b02=f22*fC20 + f23*fC30 + f24*fC40;
  Double_t b10=f02*fC21 + f03*fC31 + f04*fC41, b11=f12*fC21 + f14*fC41 + f13*fC31;
  Double_t b12=f22*fC21 + f23*fC31 + f24*fC41;
  Double_t b20=f02*fC22 + f03*fC32 + f04*fC42, b21=f12*fC22 + f14*fC42 + f13*fC32;
  Double_t b22=f22*fC22 + f23*fC32 + f24*fC42;
  Double_t b40=f02*fC42 + f03*fC43 + f04*fC44, b41=f12*fC42 + f14*fC44 + f13*fC43;
  Double_t b42=f22*fC42 + f23*fC43 + f24*fC44;
  Double_t b30=f02*fC32 + f03*fC33 + f04*fC43, b31=f12*fC32 + f14*fC43 + f13*fC33;
  Double_t b32=f22*fC32 + f23*fC33 + f24*fC43;
  
  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f04*b40+f03*b30,a01=f02*b21+f04*b41+f03*b31,a02=f02*b22+f04*b42+f03*b32;
  Double_t a11=f12*b21+f14*b41+f13*b31,a12=f12*b22+f14*b42+f13*b32;
  Double_t a22=f22*b22+f24*b42+f23*b32;

  //F*C*Ft = C + (b + bt + a)
  fC00 += b00 + b00 + a00;
  fC10 += b10 + b01 + a01; 
  fC20 += b20 + b02 + a02;
  fC30 += b30;
  fC40 += b40;
  fC11 += b11 + b11 + a11;
  fC21 += b21 + b12 + a12;
  fC31 += b31; 
  fC41 += b41;
  fC22 += b22 + b22 + a22;
  fC32 += b32;
  fC42 += b42;

  bool result=kTRUE;
  if(checkSignX&&kEndCap&&fX<0){
    result&=RotateZ(false,b);
    result&=Rotate(fAlpha-TMath::Pi(),b);
    result&=RotateZ(true,b);
  }

  return result;
}

Bool_t IlcExternalTrackParam::PropagateTo(Double_t p[3],Double_t covyz[3],Double_t covxyz[3],Double_t bz) {
  //----------------------------------------------------------------
  // Propagate this track to the plane 
  // the 3D space point "p" (with the covariance matrix "covyz" and "covxyz")
  // belongs to.
  // The magnetic field is "bz" (kG)
  //
  // The track curvature and the change of the covariance matrix
  // of the track parameters are negleted !
  // (So the "step" should be small compared with 1/curvature)
  //----------------------------------------------------------------

  Double_t f=GetSnp();
  if (TMath::Abs(f) >= kAlmost1) return kFALSE;
  Double_t r=TMath::Sqrt(1.- f*f);
  Double_t a=f/r, b=GetTgl()/r;

  Double_t s2=333.*333.;  //something reasonably big (cm^2)
 
  TMatrixDSym tV(3);
  tV(0,0)=  s2;  tV(0,1)=  a*s2;  tV(0,2)=  b*s2;
  tV(1,0)=a*s2;  tV(1,1)=a*a*s2;  tV(1,2)=a*b*s2;
  tV(2,0)=b*s2;  tV(2,1)=a*b*s2;  tV(2,2)=b*b*s2;

  TMatrixDSym pV(3);
  pV(0,0)=covxyz[0]; pV(0,1)=covxyz[1]; pV(0,2)=covxyz[2];
  pV(1,0)=covxyz[1]; pV(1,1)=covyz[0];  pV(1,2)=covyz[1];
  pV(2,0)=covxyz[2]; pV(2,1)=covyz[1];  pV(2,2)=covyz[2];

  TMatrixDSym tpV(tV);
  tpV+=pV;
  tpV.Invert();
  if (!tpV.IsValid()) return kFALSE;

  TMatrixDSym pW(3),tW(3);
  for (Int_t i=0; i<3; i++)
    for (Int_t j=0; j<3; j++) {
      pW(i,j)=tW(i,j)=0.;
      for (Int_t k=0; k<3; k++) {
	pW(i,j) += tV(i,k)*tpV(k,j);
	tW(i,j) += pV(i,k)*tpV(k,j);
      }
    }

  Double_t t[3] = {GetX(), GetY(), GetZ()};

  Double_t x=0.;
  for (Int_t i=0; i<3; i++) x += (tW(0,i)*t[i] + pW(0,i)*p[i]);  
  Double_t crv=GetC(bz);
  if (TMath::Abs(b) < kAlmost0Field) crv=0.;
  f += crv*(x-fX);
  if (TMath::Abs(f) >= kAlmost1) return kFALSE;
  fX=x;  

  fP[0]=0.;
  for (Int_t i=0; i<3; i++) fP[0] += (tW(1,i)*t[i] + pW(1,i)*p[i]);  
  fP[1]=0.;
  for (Int_t i=0; i<3; i++) fP[1] += (tW(2,i)*t[i] + pW(2,i)*p[i]);  

  return kTRUE;  
}

Double_t IlcExternalTrackParam::GetPredictedChi2(Double_t p[3],Double_t covyz[3],Double_t covxyz[3]) const {
  //----------------------------------------------------------------
  // Estimate the chi2 of the 3D space point "p" and
  // the full covariance matrix "covyz" and "covxyz"
  //
  // Cov(x,x) ... :   covxyz[0]
  // Cov(y,x) ... :   covxyz[1]  covyz[0]
  // Cov(z,x) ... :   covxyz[2]  covyz[1]  covyz[2]
  //----------------------------------------------------------------

  Double_t res[3] = {
    GetX() - p[0],
    GetY() - p[1],
    GetZ() - p[2]
  };

  Double_t f=GetSnp();
  if (TMath::Abs(f) >= kAlmost1) return kVeryBig;
  Double_t r=TMath::Sqrt(1.- f*f);
  Double_t a=f/r, b=GetTgl()/r;

  Double_t s2=333.*333.;  //something reasonably big (cm^2)
 
  TMatrixDSym v(3);
  v(0,0)=  s2;  v(0,1)=  a*s2;                 v(0,2)=  b*s2;;
  v(1,0)=a*s2;  v(1,1)=a*a*s2 + GetSigmaY2();  v(1,2)=a*b*s2 + GetSigmaZY();
  v(2,0)=b*s2;  v(2,1)=a*b*s2 + GetSigmaZY();  v(2,2)=b*b*s2 + GetSigmaZ2();

  v(0,0)+=covxyz[0]; v(0,1)+=covxyz[1]; v(0,2)+=covxyz[2];
  v(1,0)+=covxyz[1]; v(1,1)+=covyz[0];  v(1,2)+=covyz[1];
  v(2,0)+=covxyz[2]; v(2,1)+=covyz[1];  v(2,2)+=covyz[2];

  v.Invert();
  if (!v.IsValid()) return kVeryBig;

  Double_t chi2=0.;
  for (Int_t i = 0; i < 3; i++)
    for (Int_t j = 0; j < 3; j++) chi2 += res[i]*res[j]*v(i,j);

  return chi2;  


}



void IlcExternalTrackParam::Propagate(Double_t len, Double_t x[3],
Double_t p[3], Double_t bz) const {
  //+++++++++++++++++++++++++++++++++++++++++    
  // Origin: K. Shileev (Kirill.Shileev@cern.ch)
  // Extrapolate track along simple helix in magnetic field
  // Arguments: len -distance alogn helix, [cm]
  //            bz  - mag field, [kGaus]   
  // Returns: x and p contain extrapolated positon and momentum  
  // The momentum returned for straight-line tracks is meaningless !
  //+++++++++++++++++++++++++++++++++++++++++    
  GetXYZ(x);
    
  if (TMath::Abs(Get1Pt()) < kAlmost0){ //straight-line tracks
     Double_t unit[3]; GetDirection(unit);
     x[0]+=unit[0]*len;   
     x[1]+=unit[1]*len;   
     x[2]+=unit[2]*len;

     p[0]=unit[0]/kAlmost0;   
     p[1]=unit[1]/kAlmost0;   
     p[2]=unit[2]/kAlmost0;   
  } else {
     GetPxPyPz(p);
     Double_t pp=GetP();
     Double_t a = -kB2C*bz*GetSign();
     Double_t rho = a/pp;
     x[0] += p[0]*TMath::Sin(rho*len)/a - p[1]*(1-TMath::Cos(rho*len))/a;
     x[1] += p[1]*TMath::Sin(rho*len)/a + p[0]*(1-TMath::Cos(rho*len))/a;
     x[2] += p[2]*len/pp;

     Double_t p0=p[0];
     p[0] = p0  *TMath::Cos(rho*len) - p[1]*TMath::Sin(rho*len);
     p[1] = p[1]*TMath::Cos(rho*len) + p0  *TMath::Sin(rho*len);
  }
}

Bool_t IlcExternalTrackParam::PropagateToWithCorrection(Double_t xToGo, Double_t b, Double_t mass, 
Double_t maxStep, Bool_t rotateTo, Double_t maxSnp){
  //----------------------------------------------------------------
  //
  // MI's function
  //
  // Propagates this track to the plane X=xk (cm) 
  // in the magnetic field "b" (kG),
  // the correction for the material is included
  //
  // mass     - mass used in propagation - used for energy loss correction
  // maxStep  - maximal step for propagation
  //----------------------------------------------------------------
//   cout<<"Propagate to with corrections xtogo="<<xToGo<<" b "<<b<<" mass "<<mass<<" maxstep "<<
//     maxStep<<" rotateTo "<<rotateTo<<" maxsnp "<<maxSnp<<endl;
//   Print("params");
  const Double_t kEpsilon = 0.00001;
  Double_t xpos     = GetX();
  Double_t dir      = (xpos<xToGo) ? 1.:-1.;
  //
  while ( (xToGo-xpos)*dir > kEpsilon){
//     cout<<"In begin of step "<<endl;
//     Print("params");
    Double_t step = dir*TMath::Min(TMath::Abs(xToGo-xpos), maxStep);
    Double_t x    = xpos+step;
    Double_t xyz0[3],xyz1[3],param[7];
    GetXYZ(xyz0);   //starting global position
    if (!GetXYZAt(x,b,xyz1)) return kFALSE;   // no prolongation
    if(rotateTo&&((TMath::Hypot(xyz1[0],xyz1[1])-TMath::Hypot(xyz0[0],xyz0[1]))*
		  (xToGo-xpos))<0) return kFALSE;
    xyz1[2]+=kEpsilon; // waiting for bug correction in geo
    IlcKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);	
    if (TMath::Abs(GetSnpAt(x,b)) >= maxSnp) return kFALSE;
    if (!IlcExternalTrackParam::PropagateTo(x,b))  return kFALSE;
//     cout<<"After propagate "<<endl;
//     Print("params");


    Double_t rho=param[0],x0=param[1],distance=param[4];
    Double_t d=-distance*rho/x0*dir*fDir;

    if (!CorrectForMaterial(d,x0,mass)) return kFALSE;
//     cout<<"After correct "<<endl;
//     Print("params");
    if (rotateTo){
      if (TMath::Abs(GetSnp()) >= maxSnp) return kFALSE;
      GetXYZ(xyz0);   // global position
      Double_t alphan = TMath::ATan2(xyz0[1], xyz0[0]); 
      //
      double dirbefore=fDir;
      if (!Rotate(alphan,b)) return kFALSE;
      if(dirbefore!=fDir) {
	maxStep*=0.5;
	if(maxStep<=kAlmost0) return kFALSE;
      }
//       cout<<"After Rotate "<<endl;
//       Print("params");
    }
    xpos = GetX();
  }
//   cout<<"Finished"<<endl;
//   Print("params");


  return kTRUE;
}


Bool_t IlcExternalTrackParam::Intersect(Double_t pnt[3], Double_t norm[3],
Double_t bz) const {
  //+++++++++++++++++++++++++++++++++++++++++    
  // Origin: K. Shileev (Kirill.Shileev@cern.ch)
  // Finds point of intersection (if exists) of the helix with the plane. 
  // Stores result in fX and fP.   
  // Arguments: planePoint,planeNorm - the plane defined by any plane's point 
  // and vector, normal to the plane
  // Returns: kTrue if helix intersects the plane, kFALSE otherwise.
  //+++++++++++++++++++++++++++++++++++++++++    
  Double_t x0[3]; GetXYZ(x0); //get track position in MARS
  
  //estimates initial helix length up to plane
  Double_t s=
    (pnt[0]-x0[0])*norm[0] + (pnt[1]-x0[1])*norm[1] + (pnt[2]-x0[2])*norm[2];
  Double_t dist=99999,distPrev=dist;
  Double_t x[3],p[3]; 
  while(TMath::Abs(dist)>0.00001){
    //calculates helix at the distance s from x0 ALONG the helix
    Propagate(s,x,p,bz);

    //distance between current helix position and plane
    dist=(x[0]-pnt[0])*norm[0]+(x[1]-pnt[1])*norm[1]+(x[2]-pnt[2])*norm[2];

    if(TMath::Abs(dist) >= TMath::Abs(distPrev)) {return kFALSE;}
    distPrev=dist;
    s-=dist;
  }
  //on exit pnt is intersection point,norm is track vector at that point, 
  //all in MARS
  for (Int_t i=0; i<3; i++) {pnt[i]=x[i]; norm[i]=p[i];}
  return kTRUE;
}

Double_t 
IlcExternalTrackParam::GetPredictedChi2(Double_t p[2],Double_t cov[3]) const {
  //----------------------------------------------------------------
  // Estimate the chi2 of the space point "p" with the cov. matrix "cov"
  //----------------------------------------------------------------
  Double_t sdd = fC[0] + cov[0]; 
  Double_t sdz = fC[1] + cov[1];
  Double_t szz = fC[2] + cov[2];
  Double_t det = sdd*szz - sdz*sdz;

  if (TMath::Abs(det) < kAlmost0) return kVeryBig;

  Double_t d = fP[0] - p[0];
  Double_t z = (kEndCap?fX:fP[1]) - p[1];

  return (d*szz*d - 2*d*sdz*z + z*sdd*z)/det;
}

Bool_t IlcExternalTrackParam::Update(Double_t p[2], Double_t cov[3]) {
  //------------------------------------------------------------------
  // Update the track parameters with the space point "p" having
  // the covariance matrix "cov"
  //------------------------------------------------------------------
  Double_t &fP0=fP[0], &fP1=(kEndCap?fX:fP[1]), &fP2=fP[2], &fP3=fP[3], &fP4=fP[4];
  Double_t 
  &fC00=fC[0],
  &fC10=fC[1],   &fC11=fC[2],  
  &fC20=fC[3],   &fC21=fC[4],   &fC22=fC[5],
  &fC30=fC[6],   &fC31=fC[7],   &fC32=fC[8],   &fC33=fC[9],  
  &fC40=fC[10],  &fC41=fC[11],  &fC42=fC[12],  &fC43=fC[13], &fC44=fC[14];

  Double_t r00=cov[0], r01=cov[1], r11=cov[2];
  r00+=fC00; r01+=fC10; r11+=fC11;
  Double_t det=r00*r11 - r01*r01;

  if (TMath::Abs(det) < kAlmost0) return kFALSE;


  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;
 
  Double_t k00=fC00*r00+fC10*r01, k01=fC00*r01+fC10*r11;
  Double_t k10=fC10*r00+fC11*r01, k11=fC10*r01+fC11*r11;
  Double_t k20=fC20*r00+fC21*r01, k21=fC20*r01+fC21*r11;
  Double_t k30=fC30*r00+fC31*r01, k31=fC30*r01+fC31*r11;
  Double_t k40=fC40*r00+fC41*r01, k41=fC40*r01+fC41*r11;

  Double_t dy=p[0] - fP0, dz=p[1] - fP1;
  Double_t sf=fP2 + k20*dy + k21*dz;
  if (TMath::Abs(sf) > kAlmost1) {
    if(IlcDebugLevelClass()>=5) IlcDebug(5,Form("Warning in Update sin before = %f after = %f",fP2,sf));
    sf=TMath::Sign(2-TMath::Abs(sf),sf);
    if(TMath::Abs(sf) > kAlmost1) return kFALSE;
    fDir*=-1;
  }
  fP0 += k00*dy + k01*dz;
  fP1 += k10*dy + k11*dz;
  fP2  = sf;
  fP3 += k30*dy + k31*dz;
  fP4 += k40*dy + k41*dz;
  
  Double_t c01=fC10, c02=fC20, c03=fC30, c04=fC40;
  Double_t c12=fC21, c13=fC31, c14=fC41;

  fC00-=k00*fC00+k01*fC10; fC10-=k00*c01+k01*fC11;
  fC20-=k00*c02+k01*c12;   fC30-=k00*c03+k01*c13;
  fC40-=k00*c04+k01*c14; 

  fC11-=k10*c01+k11*fC11;
  fC21-=k10*c02+k11*c12;   fC31-=k10*c03+k11*c13;
  fC41-=k10*c04+k11*c14; 

  fC22-=k20*c02+k21*c12;   fC32-=k20*c03+k21*c13;
  fC42-=k20*c04+k21*c14; 

  fC33-=k30*c03+k31*c13;
  fC43-=k30*c04+k31*c14; 

  fC44-=k40*c04+k41*c14; 

  return kTRUE;
}

void 
IlcExternalTrackParam::GetHelixParameters(Double_t hlx[6], Double_t b) const {
  //--------------------------------------------------------------------
  // External track parameters -> helix parameters 
  // "b" - magnetic field (kG)
  //--------------------------------------------------------------------
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  
  hlx[0]=fP[0]; hlx[1]=fP[1]; hlx[2]=fP[2]; hlx[3]=fP[3];

  hlx[5]=fX*cs - hlx[0]*sn;               // x0
  hlx[0]=fX*sn + hlx[0]*cs;               // y0
//hlx[1]=                                 // z0
  hlx[2]=TMath::ASin(hlx[2]) + fAlpha;    // phi0
//hlx[3]=                                 // tgl
  hlx[4]=GetC(b);                         // C
}


static void Evaluate(const Double_t *h, Double_t t,
                     Double_t r[3],  //radius vector
                     Double_t g[3],  //first defivatives
                     Double_t gg[3]) //second derivatives
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives
  //--------------------------------------------------------------------
  Double_t phase=h[4]*t+h[2];
  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

  r[0] = h[5] + (sn - h[6])/h[4];
  r[1] = h[0] - (cs - h[7])/h[4];  
  r[2] = h[1] + h[3]*t;

  g[0] = cs; g[1]=sn; g[2]=h[3];
  
  gg[0]=-h[4]*sn; gg[1]=h[4]*cs; gg[2]=0.;
}

Double_t IlcExternalTrackParam::GetDCA(const IlcExternalTrackParam *p, 
Double_t b, Double_t &xthis, Double_t &xp) const {
  //------------------------------------------------------------
  // Returns the (weighed !) distance of closest approach between 
  // this track and the track "p".
  // Other returned values:
  //   xthis, xt - coordinates of tracks' reference planes at the DCA 
  //-----------------------------------------------------------
  Double_t dy2=GetSigmaY2() + p->GetSigmaY2();
  Double_t dz2=GetSigmaZ2() + p->GetSigmaZ2();
  Double_t dx2=dy2; 

  //dx2=dy2=dz2=1.;

  Double_t p1[8]; GetHelixParameters(p1,b);
  p1[6]=TMath::Sin(p1[2]); p1[7]=TMath::Cos(p1[2]);
  Double_t p2[8]; p->GetHelixParameters(p2,b);
  p2[6]=TMath::Sin(p2[2]); p2[7]=TMath::Cos(p2[2]);


  Double_t r1[3],g1[3],gg1[3]; Double_t t1=0.;
  Evaluate(p1,t1,r1,g1,gg1);
  Double_t r2[3],g2[3],gg2[3]; Double_t t2=0.;
  Evaluate(p2,t2,r2,g2,gg2);

  Double_t dx=r2[0]-r1[0], dy=r2[1]-r1[1], dz=r2[2]-r1[2];
  Double_t dm=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;

  Int_t max=27;
  while (max--) {
     Double_t gt1=-(dx*g1[0]/dx2 + dy*g1[1]/dy2 + dz*g1[2]/dz2);
     Double_t gt2=+(dx*g2[0]/dx2 + dy*g2[1]/dy2 + dz*g2[2]/dz2);
     Double_t h11=(g1[0]*g1[0] - dx*gg1[0])/dx2 + 
                  (g1[1]*g1[1] - dy*gg1[1])/dy2 +
                  (g1[2]*g1[2] - dz*gg1[2])/dz2;
     Double_t h22=(g2[0]*g2[0] + dx*gg2[0])/dx2 + 
                  (g2[1]*g2[1] + dy*gg2[1])/dy2 +
                  (g2[2]*g2[2] + dz*gg2[2])/dz2;
     Double_t h12=-(g1[0]*g2[0]/dx2 + g1[1]*g2[1]/dy2 + g1[2]*g2[2]/dz2);

     Double_t det=h11*h22-h12*h12;

     Double_t dt1,dt2;
     if (TMath::Abs(det)<1.e-33) {
        //(quasi)singular Hessian
        dt1=-gt1; dt2=-gt2;
     } else {
        dt1=-(gt1*h22 - gt2*h12)/det; 
        dt2=-(h11*gt2 - h12*gt1)/det;
     }

     if ((dt1*gt1+dt2*gt2)>0) {dt1=-dt1; dt2=-dt2;}

     //check delta(phase1) ?
     //check delta(phase2) ?

     if (TMath::Abs(dt1)/(TMath::Abs(t1)+1.e-3) < 1.e-4)
     if (TMath::Abs(dt2)/(TMath::Abs(t2)+1.e-3) < 1.e-4) {
        if ((gt1*gt1+gt2*gt2) > 1.e-4/dy2/dy2) 
	  IlcWarning(" stopped at not a stationary point !");
        Double_t lmb=h11+h22; lmb=lmb-TMath::Sqrt(fabs(lmb*lmb-4*det));
        if (lmb < 0.) 
	  IlcWarning(" stopped at not a minimum !");
        break;
     }

     Double_t dd=dm;
     for (Int_t div=1 ; ; div*=2) {
        Evaluate(p1,t1+dt1,r1,g1,gg1);
        Evaluate(p2,t2+dt2,r2,g2,gg2);
        dx=r2[0]-r1[0]; dy=r2[1]-r1[1]; dz=r2[2]-r1[2];
        dd=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
	if (dd<dm) break;
        dt1*=0.5; dt2*=0.5;
        if (div>512) {
           IlcWarning(" overshoot !"); break;
        }   
     }
     dm=dd;

     t1+=dt1;
     t2+=dt2;

  }

  if (max<=0) IlcWarning(" too many iterations !");

  Double_t cs=TMath::Cos(GetAlpha());
  Double_t sn=TMath::Sin(GetAlpha());
  xthis=r1[0]*cs + r1[1]*sn;

  cs=TMath::Cos(p->GetAlpha());
  sn=TMath::Sin(p->GetAlpha());
  xp=r2[0]*cs + r2[1]*sn;

  return TMath::Sqrt(dm*TMath::Sqrt(fabs(dy2*dz2)));
}
 
Double_t IlcExternalTrackParam::
PropagateToDCA(IlcExternalTrackParam *p, Double_t b) {
  //--------------------------------------------------------------
  // Propagates this track and the argument track to the position of the
  // distance of closest approach.
  // Returns the (weighed !) distance of closest approach.
  //--------------------------------------------------------------
  Double_t xthis,xp;
  Double_t dca=GetDCA(p,b,xthis,xp);

  if (!PropagateTo(xthis,b)) {
    //IlcWarning(" propagation failed !");
    return 1e+33;
  }

  if (!p->PropagateTo(xp,b)) {
    //IlcWarning(" propagation failed !";
    return 1e+33;
  }

  return dca;
}




Bool_t IlcExternalTrackParam::PropagateToDCA(const double *vxy, Double_t b, Double_t maxd){
  //
  // Try to relate this track to the vertex "vtx", 
  // if the (rough) transverse impact parameter is not bigger then "maxd". 
  //            Magnetic field is "b" (kG).
  //
  // a) The track gets extapolated to the DCA to the vertex.
  // b) The impact parameters and their covariance matrix are calculated.
  //
  //    In the case of success, the returned value is kTRUE
  //    (otherwise, it's kFALSE)
  //  
  Double_t alpha=GetAlpha();
  Double_t sn=TMath::Sin(alpha), cs=TMath::Cos(alpha);
  Double_t x=GetX(), y=GetParameter()[0], snp=GetParameter()[2],csp=fDir*TMath::Sqrt(1.- snp*snp);
  Double_t xv= vxy[0]*cs + vxy[1]*sn;
  Double_t yv=-vxy[0]*sn + vxy[1]*cs;
  x-=xv; y-=yv;

  //Estimate the impact parameter neglecting the track curvature
  Double_t d=TMath::Abs(x*snp - y*csp);
  if (d > maxd) return kFALSE; 

  //Propagate to the DCA
  Double_t crv=GetC(b);
  Double_t tgfv=-(crv*x - snp)/(crv*y + csp);
  sn=tgfv/TMath::Sqrt(1.+ tgfv*tgfv); cs=TMath::Sqrt(1.- sn*sn);

  x = xv*cs + yv*sn;
  yv=-xv*sn + yv*cs; xv=x;

  if (!Propagate(alpha+TMath::ASin(sn),xv,b,true)) return kFALSE;
  return kTRUE;
}




Bool_t Local2GlobalMomentum(Double_t p[3],Double_t alpha,double dir) {
  //----------------------------------------------------------------
  // This function performs local->global transformation of the
  // track momentum.
  // When called, the arguments are:
  //    p[0] = 1/pt of the track;
  //    p[1] = sine of local azim. angle of the track momentum;
  //    p[2] = tangent of the track momentum dip angle;
  //   alpha - rotation angle. 
  // The result is returned as:
  //    p[0] = px
  //    p[1] = py
  //    p[2] = pz
  // Results for (nearly) straight tracks are meaningless !
  //----------------------------------------------------------------
  if (TMath::Abs(p[0])<=kAlmost0) return kFALSE;
  if (TMath::Abs(p[1])> kAlmost1) return kFALSE;

  Double_t pt=1./TMath::Abs(p[0]);
  Double_t cs=TMath::Cos(alpha), sn=TMath::Sin(alpha);
  Double_t r=dir*TMath::Sqrt(1 - p[1]*p[1]);
  p[0]=pt*(r*cs - p[1]*sn); p[1]=pt*(p[1]*cs + r*sn); p[2]=pt*p[2];

  return kTRUE;
}

Bool_t Local2GlobalPosition(Double_t r[3],Double_t alpha) {
  //----------------------------------------------------------------
  // This function performs local->global transformation of the
  // track position.
  // When called, the arguments are:
  //    r[0] = local x
  //    r[1] = local y
  //    r[2] = local z
  //   alpha - rotation angle. 
  // The result is returned as:
  //    r[0] = global x
  //    r[1] = global y
  //    r[2] = global z
  //----------------------------------------------------------------
  Double_t cs,sn; sincos(alpha,&sn,&cs);
  Double_t x=r[0];
  r[0]=x*cs - r[1]*sn; r[1]=x*sn + r[1]*cs;

  return kTRUE;
}

void IlcExternalTrackParam::GetDirection(Double_t d[3]) const {
  //----------------------------------------------------------------
  // This function returns a unit vector along the track direction
  // in the global coordinate system.
  //----------------------------------------------------------------
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  Double_t snp=fP[2];
  Double_t csp =fDir*TMath::Sqrt(1.- snp*snp);
  Double_t norm=TMath::Sqrt(1.+ fP[3]*fP[3]);
  d[0]=(csp*cs - snp*sn)/norm; 
  d[1]=(snp*cs + csp*sn)/norm; 
  d[2]=fP[3]/norm;
}

Bool_t IlcExternalTrackParam::GetPxPyPz(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  p[0]=fP[4]; p[1]=fP[2]; p[2]=fP[3];
  return Local2GlobalMomentum(p,fAlpha,fDir);
}

Bool_t IlcExternalTrackParam::GetXYZ(Double_t *r,bool local) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  //---------------------------------------------------------------------
  r[0]=fX; r[1]=fP[0]; r[2]=fP[1];
  if(local) return kTRUE;
  return Local2GlobalPosition(r,fAlpha);
}

Bool_t IlcExternalTrackParam::GetCovarianceXYZPxPyPz(Double_t cv[21]) const {
  //---------------------------------------------------------------------
  // This function returns the global covariance matrix of the track params
  // 
  // Cov(x,x) ... :   cv[0]
  // Cov(y,x) ... :   cv[1]  cv[2]
  // Cov(z,x) ... :   cv[3]  cv[4]  cv[5]
  // Cov(px,x)... :   cv[6]  cv[7]  cv[8]  cv[9]
  // Cov(py,x)... :   cv[10] cv[11] cv[12] cv[13] cv[14]
  // Cov(pz,x)... :   cv[15] cv[16] cv[17] cv[18] cv[19] cv[20]
  //
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  if (TMath::Abs(fP[4])<=kAlmost0) {
     for (Int_t i=0; i<21; i++) cv[i]=0.;
     return kFALSE;
  }
  if (TMath::Abs(fP[2]) > kAlmost1) {
     for (Int_t i=0; i<21; i++) cv[i]=0.;
     return kFALSE;
  }
  Double_t pt=1./TMath::Abs(fP[4]);
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  Double_t r=fDir*TMath::Sqrt(1-fP[2]*fP[2]);

  Double_t m00=-sn, m10=cs;
  Double_t m23=-pt*(sn + fP[2]*cs/r);
  Double_t m43=-pt*pt*(r*cs - fP[2]*sn)*TMath::Sign(1.,fP[4]);
  Double_t m24= pt*(cs - fP[2]*sn/r);
  Double_t m44=-pt*pt*(r*sn + fP[2]*cs)*TMath::Sign(1.,fP[4]);
  Double_t m35=pt, m45=-pt*pt*fP[3]*TMath::Sign(1.,fP[4]);

  cv[0 ] = fC[0]*m00*m00;
  cv[1 ] = fC[0]*m00*m10; 
  cv[2 ] = fC[0]*m10*m10;
  cv[3 ] = fC[1]*m00; 
  cv[4 ] = fC[1]*m10; 
  cv[5 ] = fC[2];
  cv[6 ] = m00*(fC[3]*m23 + fC[10]*m43); 
  cv[7 ] = m10*(fC[3]*m23 + fC[10]*m43); 
  cv[8 ] = fC[4]*m23 + fC[11]*m43; 
  cv[9 ] = m23*(fC[5]*m23 + fC[12]*m43)  +  m43*(fC[12]*m23 + fC[14]*m43);
  cv[10] = m00*(fC[3]*m24 + fC[10]*m44); 
  cv[11] = m10*(fC[3]*m24 + fC[10]*m44); 
  cv[12] = fC[4]*m24 + fC[11]*m44; 
  cv[13] = m23*(fC[5]*m24 + fC[12]*m44)  +  m43*(fC[12]*m24 + fC[14]*m44);
  cv[14] = m24*(fC[5]*m24 + fC[12]*m44)  +  m44*(fC[12]*m24 + fC[14]*m44);
  cv[15] = m00*(fC[6]*m35 + fC[10]*m45); 
  cv[16] = m10*(fC[6]*m35 + fC[10]*m45); 
  cv[17] = fC[7]*m35 + fC[11]*m45; 
  cv[18] = m23*(fC[8]*m35 + fC[12]*m45)  +  m43*(fC[13]*m35 + fC[14]*m45);
  cv[19] = m24*(fC[8]*m35 + fC[12]*m45)  +  m44*(fC[13]*m35 + fC[14]*m45); 
  cv[20] = m35*(fC[9]*m35 + fC[13]*m45)  +  m45*(fC[13]*m35 + fC[14]*m45);

  return kTRUE;
}


Bool_t 
IlcExternalTrackParam::GetPxPyPzAt(Double_t x, Double_t b, Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum extrapolated to
  // the radial position "x" (cm) in the magnetic field "b" (kG)
  //---------------------------------------------------------------------
  p[0]=fP[4]; 
  p[1]=fP[2]+(x-fX)*GetC(b); 
  p[2]=fP[3];
  return Local2GlobalMomentum(p,fAlpha,fDir);
}

Bool_t 
IlcExternalTrackParam::GetYAt(Double_t x, Double_t b, Double_t &y) const {
  //---------------------------------------------------------------------
  // This function returns the local Y-coordinate of the intersection 
  // point between this track and the reference plane "x" (cm). 
  // Magnetic field "b" (kG)
  //---------------------------------------------------------------------
  Double_t dx=x-fX;
  if(TMath::Abs(dx)<=kAlmost0) {y=fP[0]; return kTRUE;}

  Double_t f1=fP[2], f2=f1 + dx*GetC(b);

  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;
  
  Double_t r1=fDir*TMath::Sqrt(1.- f1*f1), r2=fDir*TMath::Sqrt(1.- f2*f2);
  y = fP[0] + dx*(f1+f2)/(r1+r2);
  return kTRUE;
}

Bool_t 
IlcExternalTrackParam::GetZAt(Double_t x, Double_t b, Double_t &z) const {
  //---------------------------------------------------------------------
  // This function returns the local Z-coordinate of the intersection 
  // point between this track and the reference plane "x" (cm). 
  // Magnetic field "b" (kG)
  //---------------------------------------------------------------------
  Double_t dx=x-fX;
  if(TMath::Abs(dx)<=kAlmost0) {z=fP[1]; return kTRUE;}

  Double_t crv=GetC(b);

  Double_t f1=fP[2], f2=f1 + dx*crv;

  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;

  
  Double_t r1=sqrt(1.- f1*f1), r2=sqrt(1.- f2*f2);

  if(fabs(crv)<kAlmost0) return fP[1]+dx/r1/fDir*fP[3];
  
  double sinphi=dx*crv/2/TMath::Sqrt(0.5*(1.+r1*r2-f1*f2))/fDir;
  if(TMath::Abs(sinphi)>= kAlmost1) return kFALSE;


  z = fP[1] + 2./crv*TMath::ASin(sinphi)*fP[3];
  
  return kTRUE;
}

Bool_t IlcExternalTrackParam::GetXYZAt(Double_t x, Double_t b, Double_t *xyz,bool local) const {
  //---------------------------------------------------------------------
  // This function returns the global track position extrapolated to
  // the radial position "x" (cm) in the magnetic field "b" (kG)
  //---------------------------------------------------------------------
  if(kEndCap) return GetXYZAtZ(x,b,xyz,local);
  Double_t dx=x-fX;
  if(TMath::Abs(dx)<=kAlmost0) return GetXYZ(xyz,local);

  Double_t crv=GetC(b);

  Double_t f1=fP[2], f2=f1 + dx*crv;

  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;
  
  Double_t r1=fDir*TMath::Sqrt(1.- f1*f1), r2=fDir*TMath::Sqrt(1.- f2*f2);

  double sinphi=dx*crv/2/TMath::Sqrt(0.5*(1.+r1*r2-f1*f2))/fDir;
  if(TMath::Abs(sinphi)>= kAlmost1) return kFALSE;

  xyz[0] = x;
  xyz[1] = fP[0] + dx*(f1+f2)/(r1+r2);
  xyz[2] = fP[1] + 2./crv*TMath::ASin(sinphi)*fP[3];

  if(local) return kTRUE;
  return Local2GlobalPosition(xyz,fAlpha);
}

Bool_t IlcExternalTrackParam::GetXYZAtZ(Double_t z, Double_t b, Double_t *xyz,bool local) const {
  //---------------------------------------------------------------------
  // This function returns the global track position extrapolated to
  // the radial position "x" (cm) in the magnetic field "b" (kG)
  //---------------------------------------------------------------------
  double  dz=z-fP[1];
  if(TMath::Abs(dz)<=kAlmost0) return GetXYZ(xyz,local);
  if (TMath::Abs(fP[3]) <= kAlmost0) return kFALSE;
  Double_t f1=fP[2];
  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  Double_t r1=fDir*TMath::Sqrt(1.- f1*f1);

  Double_t crv=GetC(b);
  if(fabs(crv)<kAlmost0) crv=kAlmost0;
  Double_t dl=dz/fP[3];
  Double_t sindl,cosdl;sincos(0.5*crv*dl,&sindl,&cosdl);
  xyz[0]=fX+2./crv*sindl*(cosdl*r1-sindl*f1);
  xyz[1]=fP[0]+2./crv*sindl*(cosdl*f1+sindl*r1);
  xyz[2] =z;

  if(local) return kTRUE;
  return Local2GlobalPosition(xyz,fAlpha);
}

Bool_t IlcExternalTrackParam::GetXYZAtR(Double_t r, Double_t b, Double_t *xyz,int dir,bool local) const {
  //---------------------------------------------------------------------
  // This function returns the global track position extrapolated to
  // the radial position "r" (cm) in the magnetic field "b" (kG)
  // if dir<0 than return point in opposite direction for closest "r"
  //---------------------------------------------------------------------

  Double_t crv=GetC(b);
  if(fabs(crv)<kAlmost0Curv) crv=kAlmost0Curv;
  Double_t R=1./crv;
  Double_t f1=fP[2];
  
  if (TMath::Abs(f1) >= kAlmost1) {GetXYZ(xyz,local); return kFALSE;}

  Double_t r1=fDir*TMath::Sqrt(1.- f1*f1);

  double xcrv=fX-f1*R;
  double ycrv=fP[0]+r1*R;
  
  // direction from center of circle to origin
  double phicrv=TMath::ATan2(-ycrv,-xcrv);

  double dphi0=TMath::ATan2(fP[0]-ycrv,fX-xcrv)-phicrv;
  if(dphi0<-TMath::Pi())dphi0+=2*TMath::Pi();
  if(dphi0>TMath::Pi())dphi0-=2*TMath::Pi();
  
  double dcos=(R*R+xcrv*xcrv+ycrv*ycrv-r*r)/(2*fabs(R)*TMath::Hypot(xcrv,ycrv));

  double dsin2=(-((r*r-fX*fX-fP[0]*fP[0])*crv+2*(fX*f1-fP[0]*r1))*((r*r-fX*fX-fP[0]*fP[0])*crv+2*(fX*f1-fP[0]*r1))+
		4*r*r)/4/(xcrv*xcrv+ycrv*ycrv);

  bool isOnR=kTRUE;
  if(dsin2<0){
    isOnR=kFALSE;
    dsin2=0.;
  }
  if(dsin2>1) dsin2=1.;

  dsin2=TMath::ASin(TMath::Sqrt(dsin2));
  double dphi=(dphi0>0?1.:-1.)*(dcos>=0?dsin2:(TMath::Pi()-dsin2));

  double phi=phicrv+dphi;
  if(dir<0&&fabs(dphi)>fabs(dphi0)) phi=phicrv-dphi;
  if(dir<0&&fabs(dphi)<fabs(dphi0)) phi=phicrv+TMath::Sign(2*TMath::Pi()-fabs(dphi),dphi);

  double cphi,sphi;sincos(phi,&sphi,&cphi);
  xyz[0] =fX+fabs(R)*(cphi-f1*TMath::Sign(1.,R));
  xyz[1] =fP[0]+fabs(R)*(sphi+r1*TMath::Sign(1.,R));
  xyz[2] =fP[1]+fP[3]*(phi-phicrv-dphi0)*R;

  //  cout<<"xcrv "<<xcrv<<" "<<ycrv<<" R "<<R<<" phicrv "<<phicrv<<" dphi0 "<<dphi0<<" dphi "<<dphi<<" phi "<<phi<<" xyz "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" isOnR "<<isOnR<<" "<<dsin2<<endl;

  if(local) return isOnR;
  return isOnR&&Local2GlobalPosition(xyz,fAlpha);
}

Bool_t IlcExternalTrackParam::GetXYZAtStereo(Double_t r,Double_t tanstereoangle, Double_t b, 
					     Double_t *xyz,bool local) const {
  // not yet iplemented
  Double_t crv=GetC(b);
  if(fabs(crv)<kAlmost0Curv) crv=kAlmost0Curv;
  Double_t R=1./crv;
  Double_t f1=fP[2];
  
  if (TMath::Abs(f1) >= kAlmost1) {GetXYZ(xyz,local); return kFALSE;}

  Double_t r1=fDir*TMath::Sqrt(1.- f1*f1);

  double xcrv=fX-f1*R;
  double ycrv=fP[0]+r1*R;
  
  // direction from center of circle to origin
  double phicrv=TMath::ATan2(-ycrv,-xcrv);

  double dphi0=TMath::ATan2(fP[0]-ycrv,fX-xcrv)-phicrv;
  if(dphi0<-TMath::Pi())dphi0+=2*TMath::Pi();
  if(dphi0>TMath::Pi())dphi0-=2*TMath::Pi();
  
  double dcos=(R*R+xcrv*xcrv+ycrv*ycrv-r*r)/(2*fabs(R)*TMath::Hypot(xcrv,ycrv));

  double dsin2=(-((r*r-fX*fX-fP[0]*fP[0])*crv+2*(fX*f1-fP[0]*r1))*((r*r-fX*fX-fP[0]*fP[0])*crv+2*(fX*f1-fP[0]*r1))+
		4*r*r)/4/(xcrv*xcrv+ycrv*ycrv);

  bool isOnR=kTRUE;
  if(dsin2<0){
    isOnR=kFALSE;
    dsin2=0.;
  }
  if(dsin2>1) dsin2=1.;

  dsin2=TMath::ASin(TMath::Sqrt(dsin2));
  double dphi=(dphi0>0?1.:-1.)*(dcos>=0?dsin2:(TMath::Pi()-dsin2));

  double phi=phicrv+dphi;

  double cphi,sphi;sincos(phi,&sphi,&cphi);
  xyz[0] =fX+fabs(R)*(cphi-f1*TMath::Sign(1.,R));
  xyz[1] =fP[0]+fabs(R)*(sphi+r1*TMath::Sign(1.,R));
  xyz[2] =fP[1]+fP[3]*(phi-phicrv-dphi0)*R;
  
  if(!isOnR) {
    if(local) return kFALSE;
    Local2GlobalPosition(xyz,fAlpha);
    return kFALSE;
  }
  double alpha=tanstereoangle;
  
  double dr2=xyz[0]*xyz[0]+xyz[1]*xyz[1]-r*r-alpha*alpha*xyz[2]*xyz[2];

  double dr2dt=-2*xcrv*fabs(R)*sphi+2*ycrv*fabs(R)*cphi-2*alpha*alpha*fP[3]*R*xyz[2];
  double dr2dt2=-2*xcrv*fabs(R)*cphi-2*ycrv*fabs(R)*sphi-2*alpha*alpha*fP[3]*R*fP[3]*R;

  double D=dr2dt*dr2dt-2*dr2*dr2dt2;
  if(D<0){
    if(local) return kFALSE;
    Local2GlobalPosition(xyz,fAlpha);
    return kFALSE;
  }
  double ddphi=1e30;
  if(fabs(dr2dt2)<=kAlmost0){
    if(fabs(dr2dt2)>kAlmost0)
      ddphi=-dr2/dr2dt;
  }else{
    ddphi=(-dr2dt+TMath::Sign(TMath::Sqrt(D),dr2dt))/dr2dt2;
  }
  if(fabs(ddphi)>TMath::Pi()/4){
    if(local) return kFALSE;
    Local2GlobalPosition(xyz,fAlpha);
    return kFALSE;
  }
  
  sincos(phi+ddphi,&sphi,&cphi);
  xyz[0] =xcrv+fabs(R)*cphi;
  xyz[1] =ycrv+fabs(R)*sphi;
  xyz[2] =fP[1]+fP[3]*(phi+ddphi-phicrv-dphi0)*R;

  if(local) return kTRUE;
  return Local2GlobalPosition(xyz,fAlpha);

}

Bool_t IlcExternalTrackParam::PropagateToR(Double_t r, Double_t b){
  // propagate to radius "r" with magnetic field B(kGauss)

  double xyz[3];
  if(!GetXYZAtR(r,b,xyz)) return kFALSE;
  double angle=TMath::ATan2(xyz[1],xyz[0]);
  double x=TMath::Hypot(xyz[0],xyz[1]);
  return Propagate(angle,x,b,true);
}

Double_t IlcExternalTrackParam::GetDistance2Wire(Double_t r,Double_t phi,Double_t tanstereoangle, Double_t b,
						 double* dderiv,double* rztrack,bool localZ) const{
  //----------------------------------------------------------------
  // Get distance from track to wire "rztrack[]={wire distance, z -coordinate on wire}" 
  //                         and projection matrix dderiv
  // r-radius, phi-angle of wire on Z=0
  // stereoangle of wire
  // b - magnetic field in kGaus
  //----------------------------------------------------------------

  Double_t crv=GetC(b);
  if(fabs(crv)<kAlmost0Curv) crv=kAlmost0Curv;
  Double_t R=1./crv;
  Double_t f1=fP[2];
  
  if (TMath::Abs(f1) >= kAlmost1) {
    f1=TMath::Sign(kAlmost1,f1);
  }

  Double_t r1=fDir*TMath::Sqrt(1.- f1*f1);

  double xcrv=fX-f1*R;
  double ycrv=fP[0]+r1*R;

  double cstereo=1./TMath::Hypot(1.,tanstereoangle);
  double sstereo=cstereo*tanstereoangle;
  
  double sdangle,cdangle;sincos(phi-fAlpha,&sdangle,&cdangle);
  double r0[3]={r*cdangle,r*sdangle,0};
  double nr[3]={-sstereo*sdangle,sstereo*cdangle,cstereo};

  double t0=TMath::ATan2(f1,r1);
  double t=t0;

  double v[3]={r1,f1,fP[3]}; // *R
  double xyz[3]={fX,fP[0],fP[1]};
  double nv=nr[0]*v[0]+nr[1]*v[1]+nr[2]*v[2];

  double dt=1e30;
  int i=0;
  while(fabs(dt*R)>1e-4&&(++i<20)){
    dt=-((v[0]-nv*nr[0])*(xyz[0]-r0[0])+(v[1]-nv*nr[1])*(xyz[1]-r0[1])+
	 (v[2]-nv*nr[2])*(xyz[2]-r0[2]))/(1+v[2]*v[2]-nv*nv)*crv;

    t+=dt;
    sincos(t,&v[1],&v[0]); // *R
    xyz[0]=xcrv+v[1]*R;xyz[1]=ycrv-v[0]*R;xyz[2]=fP[1]+(t-t0)*fP[3]*R;
    nv=nr[0]*v[0]+nr[1]*v[1]+nr[2]*v[2];// *R
  }

  double nxr=nr[0]*(xyz[0]-r0[0])+nr[1]*(xyz[1]-r0[1])+nr[2]*(xyz[2]-r0[2]);

  double dl[3]={xyz[0]-r0[0]-nr[0]*nxr,xyz[1]-r0[1]-nr[1]*nxr,xyz[2]-r0[2]-nr[2]*nxr};
  
  double dln[3]={(v[1]*nr[2]-v[2]*nr[1]),(v[2]*nr[0]-v[0]*nr[2]),(v[0]*nr[1]-v[1]*nr[0])};
  double distSign=(dl[0]*dln[0]+dl[1]*dln[1]+dl[2]*dln[2])/
    TMath::Sqrt(dln[0]*dln[0]+dln[1]*dln[1]+dln[2]*dln[2]);

  distSign=fabs(distSign);

  if(dderiv){
    double deriv[5]={dl[1],dl[2],
		     -R*dl[0]-f1/r1*R*dl[1]-1./r1*fP[3]*R*dl[2],
		     (t-t0)*R*dl[2],
		     ((v[1]-f1)*dl[0]-(v[0]-r1)*dl[1]+(t-t0)*fP[3]*dl[2])*(-R)/fP[4]};
    for(int i=0;i<5;i++) {
      dderiv[i]=deriv[i];
      dderiv[i]/=distSign;
    }
  } 
  if(rztrack){
    rztrack[0]=distSign;
    rztrack[1]=nxr;
    if(!localZ) rztrack[1]*=cstereo;

    if(dderiv){
      double coeff=1./((1+v[2]*v[2]-nv*nv)+(-v[1]*dl[0]+v[0]*dl[1])*crv);
      double zder[5];
      double vec[3]={(1.+nv*nv*coeff)*nr[0]-v[0]*nv*coeff,
		     (1.+nv*nv*coeff)*nr[1]-v[1]*nv*coeff,
		     (1.+nv*nv*coeff)*nr[2]-v[2]*nv*coeff};
      zder[0]=vec[1];
      zder[1]=vec[2];
      zder[2]=vec[0]*(-R)+vec[1]*(-f1/r1*R)+vec[2]*(-fP[3]*R/r1);
      zder[3]=vec[2]*R*(t-t0)-coeff*nv*dl[2];
      zder[4]=(vec[0]*(v[1]-f1)+vec[1]*(-v[0]+r1)+vec[2]*fP[3]*(t-t0)
	       -coeff*nv*(dl[0]*v[0]+dl[1]*v[1]+dl[2]*v[2])*crv)*(-R)/fP[4];

      for(int i=0;i<5;i++) {
	dderiv[5+i]=zder[i];
	if(!localZ) dderiv[5+i]*=cstereo;
      }

//       double xphi=(xcrv-r0[0])/((xcrv-r0[0])*(xcrv-r0[0])+(ycrv-r0[1])*(ycrv-r0[1]));
//       double yphi=(ycrv-r0[1])/((xcrv-r0[0])*(xcrv-r0[0])+(ycrv-r0[1])*(ycrv-r0[1]));
      
//       double zderiv[5]=
// 	{fP[3]*R*xphi,
// 	 1,
// 	 -fP[3]*(1.+(xphi*f1-yphi*r1)*R)*R/r1,
// 	 (t-t0)*R,
// 	 fP[3]*((t-t0)+(xphi*r1+yphi*f1)*R)*(-R)/fP[4]
// 	};
//       for(int i=0;i<5;i++) {
// 	dderiv[5+i]=zderiv[i];
// 	dderiv[5+i]*=TMath::Cos(stereoangle);
//       }
    }
  }

  return distSign;
}

typedef ROOT::Math::SVector<double,5> SVector5;
typedef ROOT::Math::SVector<double,1> SVector1;
typedef ROOT::Math::SVector<double,2> SVector2;
typedef ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepSym<double,5> >  SMatrixSym5;
typedef ROOT::Math::SMatrix<double,1,1,ROOT::Math::MatRepSym<double,1> >  SMatrixSym1;
typedef ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> >  SMatrixSym2;
typedef ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepStd<double,5> >  SMatrix5;
typedef ROOT::Math::SMatrix<double,1,5,ROOT::Math::MatRepStd<double,1,5> >  SMatrix1x5;
typedef ROOT::Math::SMatrix<double,5,1,ROOT::Math::MatRepStd<double,5,1> >  SMatrix5x1;
typedef ROOT::Math::SMatrix<double,2,5,ROOT::Math::MatRepStd<double,2,5> >  SMatrix2x5;
typedef ROOT::Math::SMatrix<double,5,2,ROOT::Math::MatRepStd<double,5,2> >  SMatrix5x2;


Double_t IlcExternalTrackParam::GetPredictedChi2(Double_t r,Double_t phi,Double_t tanstereoangle, Double_t b,
						 Double_t p[2],Double_t cov[3],bool withZ,bool localZ) const{
  
  //----------------------------------------------------------------
  // Estimate the chi2 of wire "p[]={wire distance, z -coordinate on wire}" 
  //                            with the cov. matrix "cov"
  // if withZ==false use only distance to wire
  // r-radius, phi-angle of wire on Z=0
  // stereoangle of wire
  // b - magnetic field in kGaus
  //----------------------------------------------------------------

  double deriv[10],rztrack[2];
  GetDistance2Wire(r,phi,tanstereoangle,b,deriv,rztrack,localZ);
  return GetPredictedChi2(rztrack,deriv,p,cov,withZ);
}

Double_t IlcExternalTrackParam::GetPredictedChi2(Double_t rztrack[2],Double_t deriv[10],
						 Double_t p[2],Double_t cov[3],bool withZ) const{
  
  //----------------------------------------------------------------
  // Estimate the chi2 of wire "p[]={wire distance, z -coordinate on wire}" 
  //                            with the cov. matrix "cov"
  // if withZ==false use only distance to wire
  // track parameters(rztrack) on wire with projection matrix(deriv) 
  // must be calculated before by IlcExternalTrackParam::GetDistance2Wire
  //----------------------------------------------------------------

  SMatrixSym5 Covariance(fC,15);
  
  if(withZ){
    SMatrix2x5  h(deriv,10);
    SMatrixSym2 ck(cov,3);
    
    ck+=ROOT::Math::Similarity(h,Covariance);
    if(!ck.Invert()){
      Warning("Update","Problem with invertion of matrix ");
      cout<<"Track Covariance Matrix:"<<endl;
      cout<<Covariance<<endl;
      cout<<"Measurment Covariance Matrix:"<<endl;
      cout<<cov[0]<<endl;
      cout<<"Projection Matrix:"<<endl;
      cout<<h<<endl;
      return 0;
    }
    return (p[0]-rztrack[0])*(p[0]-rztrack[0])*ck(0,0)+(p[1]-rztrack[1])*(p[1]-rztrack[1])*ck(1,1)+
      2*(p[0]-rztrack[0])*(p[1]-rztrack[1])*ck(0,1);
  }else{
    SMatrix1x5  h(deriv,5);
    SMatrixSym1 ck(cov,1);
    
    ck+=ROOT::Math::Similarity(h,Covariance);
    if(!ck.Invert()){
      Warning("Update","Problem with invertion of matrix ");
      cout<<"Track Covariance Matrix:"<<endl;
      cout<<Covariance<<endl;
      cout<<"Measurment Covariance Matrix:"<<endl;
      cout<<cov[0]<<endl;
      cout<<"Projection Matrix:"<<endl;
      cout<<h<<endl;
      return 0;
    }
    return (p[0]-rztrack[0])*(p[0]-rztrack[0])*ck(0,0);
  }
}

Bool_t IlcExternalTrackParam::UpdateWithWire(Double_t r,Double_t phi,Double_t tanstereoangle, Double_t b,
					     Double_t p[2], Double_t cov[3],bool withZ,bool localZ) {
  
  //----------------------------------------------------------------
  // Update the track parameters with wire mesuarement "p[]={wire distance, z -coordinate on wire}" 
  //                                                with the cov. matrix "cov"
  // if withZ==false use only distance to wire
  // r-radius, phi-angle of wire on Z=0
  // stereoangle of wire
  // b - magnetic field in kGauss
  //----------------------------------------------------------------

  double deriv[10],rztrack[2];
  GetDistance2Wire(r,phi,tanstereoangle,b,deriv,rztrack,localZ);
  return UpdateWithWire(rztrack,deriv,p,cov,withZ);
}

Bool_t IlcExternalTrackParam::UpdateWithWire(Double_t rztrack[2], Double_t deriv[10],
					     Double_t p[2], Double_t cov[3],bool withZ) {
  //----------------------------------------------------------------
  // Update the track parameters with wire mesuarement "p[]={wire distance, z -coordinate on wire}" 
  //                                        with the cov. matrix "cov"
  // if withZ==false use only distance to wire
  // track parameters(rztrack) on wire with projection matrix(deriv) 
  // must be calculated before by IlcExternalTrackParam::GetDistance2Wire
  //----------------------------------------------------------------

  SMatrixSym5 Covariance(fC,15);
  SVector5 dstate;

  if(withZ){
    SMatrix2x5  h(deriv,10);
    SMatrixSym2 ck(cov,3);
    
    ck+=ROOT::Math::Similarity(h,Covariance);
    
    if(!ck.Invert()){
      Warning("Update","Problem with invertion of matrix ");
      cout<<"Track Covariance Matrix:"<<endl;
      cout<<Covariance<<endl;
      cout<<"Measurment Covariance Matrix:"<<endl;
      cout<<cov[0]<<endl;
      cout<<"Projection Matrix:"<<endl;
      cout<<h<<endl;
      return 0;
    }
    SMatrix5x2 k=Covariance*ROOT::Math::Transpose(h)*ck;
    Covariance-=ROOT::Math::Similarity(Covariance,ROOT::Math::SimilarityT(h,ck));
    SVector2 m;m[0]=p[0]-rztrack[0];m[1]=p[1]-rztrack[1];
    dstate=k*m;
  }else{
    SMatrix1x5  h(deriv,5);
    SMatrixSym1 ck(cov,1);
    
    ck+=ROOT::Math::Similarity(h,Covariance);
    
    if(!ck.Invert()){
      Warning("Update","Problem with invertion of matrix ");
      cout<<"Track Covariance Matrix:"<<endl;
      cout<<Covariance<<endl;
      cout<<"Measurment Covariance Matrix:"<<endl;
      cout<<cov[0]<<endl;
      cout<<"Projection Matrix:"<<endl;
      cout<<h<<endl;
      return 0;
    }
    SMatrix5x1 k=Covariance*ROOT::Math::Transpose(h)*ck;
    Covariance-=ROOT::Math::Similarity(Covariance,ROOT::Math::SimilarityT(h,ck));
    SVector1 m;m[0]=p[0]-rztrack[0];
    dstate=k*m;
  }

  //  cout<<"mes "<<rztrack<<" "<<p[0]<<" "<<1./ck(0,0)<<" "<<cov[0]<<endl;
  for(int i=0;i<5;i++) fP[i]+=dstate[i];
  if(fabs(fP[2])>kAlmost1) {
    int nover=int(fabs(fP[2])+1)/2;
    fP[2]-=TMath::Sign(nover*2.,fP[2]);
    if(nover%2==1){
      fP[2]*=-1;
      fDir*=-1;
    }
    fP[2]=TMath::Min(fabs(fP[2]),kAlmost1)*TMath::Sign(1.,fP[2]);
  }
  for(int i=0;i<15;i++) fC[i]=Covariance.Array()[i]; 
  return 1;
}

//_____________________________________________________________________________
void IlcExternalTrackParam::Print(Option_t* option) const
{
// print the parameters and the covariance matrix

  printf("IlcExternalTrackParam: x = %-12g  alpha = %-12g\n", fX, fAlpha);
  printf("  parameters: %12g %12g %12g %12g %12g dir=%3g isendcap=%i\n",
	 fP[0], fP[1], fP[2], fP[3], fP[4],fDir,int(IsEndCap()));
  TString opt(option);
  if(!opt.Contains("param")){
    printf("  covariance: %12g\n", fC[0]);
    printf("              %12g %12g\n", fC[1], fC[2]);
    printf("              %12g %12g %12g\n", fC[3], fC[4], fC[5]);
    printf("              %12g %12g %12g %12g\n", 
	   fC[6], fC[7], fC[8], fC[9]);
    printf("              %12g %12g %12g %12g %12g\n", 
	   fC[10], fC[11], fC[12], fC[13], fC[14]);
  }
}

Double_t IlcExternalTrackParam::GetSnpAt(Double_t x,Double_t b) const {
  //
  // Get sinus at given x
  //
  Double_t crv=GetC(b);
  if (TMath::Abs(b) < kAlmost0Field) crv=0.;
  Double_t dx = x-fX;
  Double_t res = fP[2]+dx*crv;
  return res;
}
