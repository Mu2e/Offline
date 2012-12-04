#ifndef ILCEXTERNALTRACKPARAM_H
#define ILCEXTERNALTRACKPARAM_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

/* $Id: IlcExternalTrackParam.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

/*****************************************************************************
 *              "External" track parametrisation class                       *
 *                                                                           *
 *      external param0:   local Y-coordinate of a track (cm)                *
 *      external param1:   local Z-coordinate of a track (cm)                *
 *      external param2:   local sine of the track momentum azimuthal angle  *
 *      external param3:   tangent of the track momentum dip angle           *
 *      external param4:   1/pt (1/(GeV/c))                                  *
 *                                                                           *
 * The parameters are estimated at an exact position x in a local coord.     *
 * system rotated by angle alpha with respect to the global coord.system.    *
 *        Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                     *
 *****************************************************************************/
#include "TObject.h"
#include "TMath.h"

const Double_t kAlmost1=0.99999999;
const Double_t kAlmost0=1e-33;
const Double_t kVeryBig=1./kAlmost0;
const Double_t kSmallPt=0.001;

const Double_t kB2C=0.299792458e-3;
const Double_t kAlmost0Curv = 2*1e-4/1e6;
const Double_t kAlmost0Field=1.e-9;
const Double_t kMostProbablePt=0.35;

//class IlcESDVertex;

inline double Phi_mpi_pi(double phi){return phi-TMath::Nint(phi/2*TMath::InvPi())*TMath::TwoPi();}
inline double Phi_0_2pi(double phi){return phi-TMath::Nint(phi/2*TMath::InvPi()-0.5)*TMath::TwoPi();}

Double_t ApproximateBetheBloch(Double_t);

class IlcExternalTrackParam: public TObject {
 public:
  IlcExternalTrackParam();
  IlcExternalTrackParam(const IlcExternalTrackParam &);
  IlcExternalTrackParam(Double_t x, Double_t alpha, 
			const Double_t param[5], const Double_t covar[15],
			Double_t fDirection=1);
  virtual ~IlcExternalTrackParam(){}

  void Set(Double_t x,Double_t alpha,
	   const Double_t param[5], const Double_t covar[15],Double_t fDirection=1.);
  void Set(const IlcExternalTrackParam *t);
  void Reset();
  void ResetCovariance(Double_t s2,bool savecorrelation=true) {
    double cr=(savecorrelation?s2:0.);
    fC[0] *= s2;
    fC[1] *= cr;  fC[2] *= s2;
    fC[3] *= cr;  fC[4] *= cr;  fC[5] *= s2;
    fC[6] *= cr;  fC[7] *= cr;  fC[8] *= cr;  fC[9] *= s2;
    fC[10]*= cr;  fC[11]*= cr;  fC[12]*= cr;  fC[13]*= cr;  fC[14]*=s2;
  }

  const Double_t *GetParameter() const {return fP;}
  const Double_t *GetCovariance() const {return fC;}

  Double_t GetAlpha() const {return fAlpha;}
  Double_t GetX() const {return kEndCap?fP[1]:fX;}
  Double_t GetY()    const {return fP[0];}
  Double_t GetZ()    const {return kEndCap?fX:fP[1];}
  Double_t GetSnp()  const {return fP[2];}
  Double_t GetTgl()  const {return fP[3];}
  Double_t Get1Pt()  const {return fP[4];}

  Double_t GetSigmaY2() const {return fabs(fC[0]);}
  Double_t GetSigmaZY() const {return fC[1];}
  Double_t GetSigmaZ2() const {return fabs(fC[2]);}
  Double_t GetSigmaSnpY() const {return fC[3];}
  Double_t GetSigmaSnpZ() const {return fC[4];}
  Double_t GetSigmaSnp2() const {return fabs(fC[5]);}
  Double_t GetSigmaTglY() const {return fC[6];}
  Double_t GetSigmaTglZ() const {return fC[7];}
  Double_t GetSigmaTglSnp() const {return fC[8];}
  Double_t GetSigmaTgl2() const {return fabs(fC[9]);}
  Double_t GetSigma1PtY() const {return fC[10];}
  Double_t GetSigma1PtZ() const {return fC[11];}
  Double_t GetSigma1PtSnp() const {return fC[12];}
  Double_t GetSigma1PtTgl() const {return fC[13];}
  Double_t GetSigma1Pt2() const {return fabs(fC[14]);}

  Double_t GetSign() const {return (fP[4]>0) ? 1 : -1;}
  Double_t GetP() const;
  Double_t GetPt() const {
    return (TMath::Abs(fP[4])>kAlmost0) ? 1./fP[4]:TMath::Sign(kVeryBig,fP[4]);
  }
  Double_t Get1P() const;
  Double_t GetC(Double_t b) const {return fP[4]*b*kB2C;}
  void GetDZ(Double_t x,Double_t y,Double_t z,Double_t b,Float_t dz[2]) const; 
  Double_t GetD(Double_t xv, Double_t yv, Double_t b) const; 
  Double_t GetLinearD(Double_t xv, Double_t yv) const; 
  Bool_t CorrectForMaterial(Double_t d, Double_t x0, Double_t mass,
			    Double_t (*f)(Double_t)=ApproximateBetheBloch);
  Double_t GetPredictedChi2(Double_t p[2],Double_t cov[3]) const;
  Bool_t Update(Double_t p[2],Double_t cov[3]);

  // forr wires 
  Bool_t UpdateWithWire(Double_t rztrack[2], Double_t deriv[10],
			Double_t p[2], Double_t cov[3],bool withZ=false);
  Bool_t UpdateWithWire(Double_t r,Double_t phi,Double_t tanstereoangle, Double_t b,
			Double_t p[2], Double_t cov[3],bool withZ=false,bool localZ=true);
  Double_t GetPredictedChi2(Double_t rztrack[2],Double_t deriv[10],
			    Double_t p[2],Double_t cov[3],bool withZ=false) const;
  Double_t GetPredictedChi2(Double_t r,Double_t phi,Double_t tanstereoangle, Double_t b,
			    Double_t p[2],Double_t cov[3],bool withZ=false,bool localZ=true) const;
  Double_t GetDistance2Wire(Double_t r,Double_t phi,Double_t tanstereoangle, Double_t b,
			    Double_t* dderiv=0,Double_t* rztrack=0,bool localZ=true) const;
  Bool_t PropagateToR(Double_t r, Double_t b);
  Bool_t GetXYZAtStereo(Double_t r, Double_t tanstereoangle,Double_t b, Double_t xyz[3],bool local=kFALSE) const;

  bool ChangeDirection();

  Bool_t Rotate(Double_t alpha,Double_t bz);
  Bool_t RotateZ(bool toendcap,Double_t bz);
  Bool_t PropagateTo(Double_t x, Double_t b,bool byR=false,bool checkSignX=true);
  Bool_t Propagate(Double_t alpha, Double_t x, Double_t b,bool byR=false) {
    if (Rotate(alpha,b))
      if (PropagateTo(x,b,byR)) return kTRUE;
    return kFALSE;
  }
  Double_t GetNLoopsToZ(Double_t z,Double_t b);

  Double_t GetPredictedChi2(Double_t p[3],Double_t covyz[3],Double_t covxyz[3]) const;
  Bool_t PropagateTo(Double_t p[3],Double_t covyz[3],Double_t covxyz[3],Double_t b);

  void   Propagate(Double_t len,Double_t x[3],Double_t p[3],Double_t bz) const;

  Bool_t PropagateToWithCorrection(Double_t xToGo, Double_t b, Double_t mass, 
				   Double_t maxStep, Bool_t rotateTo, Double_t maxSnp);

  Bool_t Intersect(Double_t pnt[3], Double_t norm[3], Double_t bz) const;

  void GetHelixParameters(Double_t h[6], Double_t b) const;
  Double_t GetDCA(const IlcExternalTrackParam *p, Double_t b,
		  Double_t &xthis,Double_t &xp) const;
  Double_t PropagateToDCA(IlcExternalTrackParam *p, Double_t b);
  Bool_t PropagateToDCA(const double *vxy, Double_t b, Double_t maxd=1e30);

  void GetDirection(Double_t d[3]) const;
  Bool_t GetPxPyPz(Double_t *p) const;
  Double_t GetPhi() const {return fAlpha+TMath::ATan2(fP[0],fX);};
  Double_t GetR() const {return TMath::Hypot(fP[0],fX);};
  Bool_t GetXYZ(Double_t *p,bool local=kFALSE) const;
  Bool_t GetCovarianceXYZPxPyPz(Double_t cv[21]) const;
  Bool_t GetPxPyPzAt(Double_t x, Double_t b, Double_t p[3]) const;
  Bool_t GetXYZAt(Double_t x, Double_t b, Double_t xyz[3],bool local=kFALSE) const;
  Bool_t GetXYZAtZ(Double_t z, Double_t b, Double_t xyz[3],bool local=kFALSE) const;
  Bool_t GetXYZAtR(Double_t r, Double_t b, Double_t xyz[3],int dir=0,bool local=kFALSE) const;
  Bool_t GetYAt(Double_t x,  Double_t b,  Double_t &y) const;
  Bool_t GetZAt(Double_t x,  Double_t b,  Double_t &z) const;
  void Print(Option_t* option = "") const;
  Double_t GetSnpAt(Double_t x,Double_t b) const;

  void SetEndCap(bool flag=true){kEndCap=flag;}
  bool IsEndCap() const {return kEndCap;}
  Double_t GetDir() const {return fDir;}
  Double_t GetCosRDir() const {return fDir*TMath::Sqrt(1-fP[2]*fP[2])*fX+fP[2]*fP[0];}

  Double_t &Par(Int_t i) {return fP[i];}
  Double_t &Cov(Int_t i) {return fC[i];}  
  Double_t &Cov(int i,int j) {return fC[j<=i?(i*(i+1)/2+j):(j*(j+1)/2+i)];}
private:
  Bool_t               kEndCap; // Propagation in endcap in Z -direction
  Double_t             fX;     // X coordinate for the point of parametrisation
  Double_t             fAlpha; // Local <-->global coor.system rotation angle
  Double_t             fP[5];  // The track parameters
  Double_t             fDir;   // fP[5]- sign of cos(phi)
  Double_t             fC[15]; // The track parameter covariance matrix

  ClassDef(IlcExternalTrackParam, 7)
};

#endif
