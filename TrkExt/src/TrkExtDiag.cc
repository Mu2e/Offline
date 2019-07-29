//
// $Id: TrkExtDiag.cc,v 1.4 2013/05/16 18:23:39 mjlee Exp $
// $Author: mjlee $ 
// $Date: 2013/05/16 18:23:39 $
//
// Functions for reading TrkExt
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/fwd.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

#include "BFieldGeom/inc/BFieldManager.hh"
#include "TrkExt/inc/TrkExtDiag.hh"
#include "TrkExt/inc/TrkExtMCHits.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/TrkHit.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/BField/BField.hh"
#include "BTrk/BField/BFieldFixed.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/BbrPointErr.hh"
#include "BTrk/BbrGeom/BbrError.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "RecoDataProducts/inc/TrkExtTraj.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/PointTrajectory.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "TTree.h"

using namespace std; 
using CLHEP::HepVector; 
using CLHEP::Hep3Vector; 
 
namespace mu2e 
{

  TrkExtDiag::TrkExtDiag(fhicl::ParameterSet const& pset) :
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel","makeSH")),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel","g4run"))
  { }

  TrkExtDiag::~TrkExtDiag() {
  }

  TTree* TrkExtDiag::createTrkExtDiag() {
    art::ServiceHandle<art::TFileService> tfs;

    _extdiag = tfs->make<TTree>("hTrk", "Track and recon info");
    //event-wide
    _extdiag->Branch ("exitcode", &_exitcode, "exitcode/I");
    //hit
    _extdiag->Branch ("nhots", &_nhots, "nhots/i");
    _extdiag->Branch ("hotx", _hotx, "hotx[nhots]/F");
    _extdiag->Branch ("hoty", _hoty, "hoty[nhots]/F");
    _extdiag->Branch ("hotz", _hotz, "hotz[nhots]/F");
    _extdiag->Branch ("hott0", _hott0, "hott0[nhots]/F");
    //trk
    _extdiag->Branch ("ntrks", &_ntrks, "ntrks/i");
    _extdiag->Branch ("trkl0", &_trkl0, "trkl0/F");
    _extdiag->Branch ("trkl1", &_trkl1, "trkl1/F");
    _extdiag->Branch ("trkl", _trkl, "trkl[ntrks]/F");
    _extdiag->Branch ("trkx", _trkx, "trkx[ntrks]/F");
    _extdiag->Branch ("trky", _trky, "trky[ntrks]/F");
    _extdiag->Branch ("trkz", _trkz, "trkz[ntrks]/F");
    _extdiag->Branch ("trkpx", _trkpx, "trkpx[ntrks]/F");
    _extdiag->Branch ("trkpy", _trkpy, "trkpy[ntrks]/F");
    _extdiag->Branch ("trkpz", _trkpz, "trkpz[ntrks]/F");
    _extdiag->Branch ("trkp", _trkp, "trkp[ntrks]/F");
    _extdiag->Branch ("trk0", &_trk0, "d0/F:p0/F:om/F:z0/F:td/F:ed0/F:ep0/F:eom/F:ez0/F:etd/F");
    _extdiag->Branch ("trk1", &_trk1, "d0/F:p0/F:om/F:z0/F:td/F:ed0/F:ep0/F:eom/F:ez0/F:etd/F");
    //mc
    _extdiag->Branch ("nsim", &_nsim, "nsim/i");
    _extdiag->Branch ("simx", _simx, "simx[nsim]/F");
    _extdiag->Branch ("simy", _simy, "simy[nsim]/F");
    _extdiag->Branch ("simz", _simz, "simz[nsim]/F");
    _extdiag->Branch ("simp0", &_simp0, "simp0/F");
    _extdiag->Branch ("simt0", &_simt0, "simt0/F");
    //mc - turn around point
    _extdiag->Branch ("simtp", &_simtp, "x/F:y:z");
    _extdiag->Branch ("simtpqual", &_simtpqual, "simtpqual/I");
    //mc - pa and st
    _extdiag->Branch ("nmcpa", &_nmcpa, "nmcpa/i");
    _extdiag->Branch ("nmcst", &_nmcst, "nmcst/i");
    _extdiag->Branch ("mcpapx", _mcpapx, "mcpapx[nmcpa]/F");
    _extdiag->Branch ("mcpapy", _mcpapy, "mcpapy[nmcpa]/F");
    _extdiag->Branch ("mcpapz", _mcpapz, "mcpapz[nmcpa]/F");
    _extdiag->Branch ("mcpap", _mcpap, "mcpap[nmcpa]/F");
    _extdiag->Branch ("mcpadp", _mcpadp, "mcpadp[nmcpa]/F");
    _extdiag->Branch ("mcpade", _mcpade, "mcpade[nmcpa]/F");
    _extdiag->Branch ("mcpadei", _mcpadei, "mcpadei[nmcpa]/F");
    _extdiag->Branch ("mcpadeni", _mcpadeni, "mcpadeni[nmcpa]/F");
    _extdiag->Branch ("mcpaz", _mcpaz, "mcpaz[nmcpa]/F");
    _extdiag->Branch ("mcstpx", _mcstpx, "mcstpx[nmcst]/F");
    _extdiag->Branch ("mcstpy", _mcstpy, "mcstpy[nmcst]/F");
    _extdiag->Branch ("mcstpz", _mcstpz, "mcstpz[nmcst]/F");
    _extdiag->Branch ("mcstp", _mcstp, "mcstp[nmcst]/F");
    _extdiag->Branch ("mcstdp", _mcstdp, "mcstdp[nmcst]/F");
    _extdiag->Branch ("mcstde", _mcstde, "mcstde[nmcst]/F");
    _extdiag->Branch ("mcstdei", _mcstdei, "mcstdei[nmcst]/F");
    _extdiag->Branch ("mcstdeni", _mcstdeni, "mcstdeni[nmcst]/F");
    _extdiag->Branch ("mcstx", _mcstx, "mcstx[nmcst]/F");
    _extdiag->Branch ("mcsty", _mcsty, "mcsty[nmcst]/F");
    _extdiag->Branch ("mcstz", _mcstz, "mcstz[nmcst]/F");
    _extdiag->Branch ("mcstt", _mcstt, "mcstt[nmcst]/F");
    // mc - vd
    _extdiag->Branch ("vdsi", &_vdsi, "x[2]/F:y[2]/F:z[2]/F:px[2]/F:py[2]/F:pz[2]/F:p[2]/F:d0[2]/F:p0[2]/F:om[2]/F:z0[2]/F:td[2]/F:status/i");
    _extdiag->Branch ("vdso", &_vdso, "x[2]/F:y[2]/F:z[2]/F:px[2]/F:py[2]/F:pz[2]/F:p[2]/F:d0[2]/F:p0[2]/F:om[2]/F:z0[2]/F:td[2]/F:status/i");
    _extdiag->Branch ("vdtf", &_vdtf, "x[2]/F:y[2]/F:z[2]/F:px[2]/F:py[2]/F:pz[2]/F:p[2]/F:d0[2]/F:p0[2]/F:om[2]/F:z0[2]/F:td[2]/F:status/i");
    _extdiag->Branch ("vdtm", &_vdtm, "x[2]/F:y[2]/F:z[2]/F:px[2]/F:py[2]/F:pz[2]/F:p[2]/F:d0[2]/F:p0[2]/F:om[2]/F:z0[2]/F:td[2]/F:status/i");
    _extdiag->Branch ("vdtb", &_vdtb, "x[2]/F:y[2]/F:z[2]/F:px[2]/F:py[2]/F:pz[2]/F:p[2]/F:d0[2]/F:p0[2]/F:om[2]/F:z0[2]/F:td[2]/F:status/i");
    // trkext
    _extdiag->Branch ("next", &_next, "next/i");
    _extdiag->Branch ("extx", _extx, "extx[next]/F");
    _extdiag->Branch ("exty", _exty, "exty[next]/F");
    _extdiag->Branch ("extz", _extz, "extz[next]/F");
    _extdiag->Branch ("extpx", _extpx, "extpx[next]/F");
    _extdiag->Branch ("extpy", _extpy, "extpy[next]/F");
    _extdiag->Branch ("extpz", _extpz, "extpz[next]/F");
    _extdiag->Branch ("extp", _extp, "extp[next]/F");
    _extdiag->Branch ("extrho", _extrho, "extrho[next]/F");
    _extdiag->Branch ("exts", _exts, "exts[next]/F");
    _extdiag->Branch ("extt", _extt, "extt[next]/F");
    _extdiag->Branch ("extex", _extex, "extex[next]/F");
    _extdiag->Branch ("extey", _extey, "extey[next]/F");
    _extdiag->Branch ("extez", _extez, "extez[next]/F");
    _extdiag->Branch ("extepx", _extepx, "extepx[next]/F");
    _extdiag->Branch ("extepy", _extepy, "extepy[next]/F");
    _extdiag->Branch ("extepz", _extepz, "extepz[next]/F");
    _extdiag->Branch ("extep", _extep, "extep[next]/F");
    _extdiag->Branch ("exter", _exter, "exter[next]/F");
    _extdiag->Branch ("extvid", _extvid, "extvid[next]/I");
    // trkext - pa and st
    _extdiag->Branch ("nextpa", &_nextpa, "nextpa/i");
    _extdiag->Branch ("nextst", &_nextst, "nextst/i");
    _extdiag->Branch ("extpaz", _extpaz, "extpaz[nextpa]/F");
    _extdiag->Branch ("extpadp",&_extpadp, "extpadp[nextpa]/F");
    _extdiag->Branch ("extpadptot", &_extpadptot, "extpadptot/F");
    _extdiag->Branch ("extstx", _extstx, "extstx[nextst]/F");
    _extdiag->Branch ("extsty", _extsty, "extsty[nextst]/F");
    _extdiag->Branch ("extstz", _extstz, "extstz[nextst]/F");
    _extdiag->Branch ("extstt", _extstt, "extstt[nextst]/F");
    _extdiag->Branch ("extstdp", _extstdp, "extstdp[nextst]/F");
    _extdiag->Branch ("extstdptot", &_extstdptot, "extstdptot/F");
    // trkext - vd
    _extdiag->Branch ("extvdsi", &_extvdsi, "r[2]/F:p[2]/F:d0[2]/F:p0[2]/F:om[2]/F:z0[2]/F:td[2]/F:er[2]/F:ep[2]/F:ed0[2]/F:ep0[2]/F:eom[2]/F:ez0[2]/F:etd[2]/F:status/i"); 
    _extdiag->Branch ("extvdso", &_extvdso, "r[2]/F:p[2]/F:d0[2]/F:p0[2]/F:om[2]/F:z0[2]/F:td[2]/F:er[2]/F:ep[2]/F:ed0[2]/F:ep0[2]/F:eom[2]/F:ez0[2]/F:etd[2]/F:status/i"); 
    _extdiag->Branch ("extvdtf", &_extvdtf, "r[2]/F:p[2]/F:d0[2]/F:p0[2]/F:om[2]/F:z0[2]/F:td[2]/F:er[2]/F:ep[2]/F:ed0[2]/F:ep0[2]/F:eom[2]/F:ez0[2]/F:etd[2]/F:status/i"); 
    _extdiag->Branch ("extvdtm", &_extvdtm, "r[2]/F:p[2]/F:d0[2]/F:p0[2]/F:om[2]/F:z0[2]/F:td[2]/F:er[2]/F:ep[2]/F:ed0[2]/F:ep0[2]/F:eom[2]/F:ez0[2]/F:etd[2]/F:status/i"); 
    _extdiag->Branch ("extvdtb", &_extvdtb, "r[2]/F:p[2]/F:d0[2]/F:p0[2]/F:om[2]/F:z0[2]/F:td[2]/F:er[2]/F:ep[2]/F:ed0[2]/F:ep0[2]/F:eom[2]/F:ez0[2]/F:etd[2]/F:status/i"); 
    // trkext - turn around point
    _extdiag->Branch ("extitp", &_extitp, "extitp/I");
    _extdiag->Branch ("exttp", &_exttp, "x/F:y:z:px:py:pz:t:s:volid/I");
;
    return _extdiag;
  }


  void TrkExtDiag::setRunInfo(){ 
    GeomHandle<DetectorSystem> det;
    _origin = det->toMu2e( CLHEP::Hep3Vector(0.,0.,0.) );   // add this to transfer detector coord. to mu2e coord.
    _vdzsi = _vdzso = _vdztf = _vdztm = _vdztb = 0;
    GeomHandle<VirtualDetector> vd;
    _vdzsi = (vd->getGlobal(VirtualDetectorId::ST_In) - _origin).z();
    _vdzso = (vd->getGlobal(VirtualDetectorId::ST_Out) - _origin).z();
    _vdztf = (vd->getGlobal(VirtualDetectorId::TT_FrontHollow) - _origin).z();
    _vdztm = (vd->getGlobal(VirtualDetectorId::TT_MidInner) - _origin).z();
    _vdztb = (vd->getGlobal(VirtualDetectorId::TT_Back) - _origin).z();
  }

  void TrkExtDiag::setSubRunInfo(){ 
    _bfMgr = GeomHandle<BFieldManager>().get();
  }

  HepVector TrkExtDiag::getHelixParameters (const Hep3Vector & x, const Hep3Vector & p, int sign) const {
    HepVector hpar(6);
    const HepPoint vertex(x.x(), x.y(), x.z());
    double Bz = (_bfMgr->getBField(x + _origin)).z();
    double fltlen = 0.;
    TrkHelixUtils::helixFromMom(hpar, fltlen, vertex, p, -1.*(double)sign, Bz);
    return hpar;
  }

  HepVector TrkExtDiag::getHelixParametersErr (const Hep3Vector & x, const Hep3Vector & p, Hep3Vector & ex, Hep3Vector & ep, int sign) const {
    HepVector hpar(6);
    const HepPoint vertex(x.x(), x.y(), x.z());
    Hep3Vector B = _bfMgr->getBField(x + _origin);
    //double fltlen = 0.;
    BbrError exm(3);
    BbrError epm(3);
    exm(1,1) = ex.x();
    exm(2,2) = ex.y();
    exm(3,3) = ex.z();
    epm(1,1) = ep.x();
    epm(2,2) = ep.y();
    epm(3,3) = ep.z();
    BbrPointErr pos(vertex, exm);
    BbrVectorErr pmom(p, epm);
    HepMatrix cxp(3,3,0);
    BFieldFixed fieldmap (B.x(), B.y(), B.z());

    HelixParams hepar = TrkHelixUtils::helixFromMomErr ( pos, pmom, cxp, -1.*(double) sign, fieldmap) ;
    const HepSymMatrix & cov = hepar.covariance();
    for (int i = 1 ; i <= 6 ; ++i) {
      hpar(i) = sqrt(cov(i,i));
    }
    return hpar;
  }

  double TrkExtDiag::findTurnAroundSim (int i, double f1, double f2, double f3) {
    double a = 0.5*f1 - f2 + 0.5*f3;
    double b = -0.5*f1*(2*i+1) + f2*(2*i) - 0.5*f3*(2*i-1);
    if (a!=0) return -b/2./a; 
    else return 0;
  }

  double TrkExtDiag::findTurnAround (double s1, double s2, double s3, double f1, double f2, double f3) {
    double s1s2 = s1-s2;
    double s2s3 = s2-s3;
    double s3s1 = s3-s1;
    double a = -f1/s1s2/s3s1 - f2/s1s2/s2s3 - f3/s3s1/s2s3;
    double b = +f1*(s2+s3)/s1s2/s3s1 + f2*(s1+s3)/s1s2/s2s3 + f3*(s1+s2)/s3s1/s2s3;
    double c = -f1*s2*s3/s1s2/s3s1 - f2*s1*s3/s1s2/s2s3 - f3*s1*s2/s3s1/s2s3;
    double D = b*b-4.*a*c;
    if (D>=0) {
      double r1 = (-b+sqrt(D))/2./a;
      double r2 = (-b-sqrt(D))/2./a;
      double d1 = fabs(r1-s2);
      double d2 = fabs(r2-s2);
      if (d1>d2) return r2;
      else if (d2>d1) return r1;
      else return 0.5*(r1+r2);
    }
    else {
      return 0; 
    }
  }

  double TrkExtDiag::interpolate2 (double s, double s1, double s2, double f1, double f2) {
    return (f2-f1)/(s2-s1)*(s-s1)+f1;
  }

  double TrkExtDiag::interpolate (double s, double s1, double s2, double s3, double f1, double f2, double f3) {
    double s1s2 = s1-s2;
    double s2s3 = s2-s3;
    double s3s1 = s3-s1;
    return -f1*(s-s2)*(s-s3)/s1s2/s3s1 - f2*(s-s1)*(s-s3)/s1s2/s2s3 - f3*(s-s1)*(s-s2)/s3s1/s2s3;
  }

  double TrkExtDiag::getRadialError (const CLHEP::Hep3Vector & xx, const CLHEP::HepMatrix & cov) {
    double r = xx.rho();
    double x = xx.x();
    double y = xx.y();
    if (r == 0) r = 0.001;
    //cout << xx<<endl;
    //cout << cov<<endl;
    return sqrt(x*x*cov(1,1)+y*y*cov(2,2)+2.*x*y*cov(1,2))/r;
  }

  double TrkExtDiag::getMomentumError (const CLHEP::Hep3Vector & p, const CLHEP::HepMatrix & cov) {
    double pr = p.mag();
    double px = p.x();
    double py = p.y();
    double pz = p.z();
    if (pr == 0) pr = 0.001;
    return sqrt(px*px*cov(4,4)+py*py*cov(5,5)+pz*pz*cov(6,6)+2.*px*py*cov(4,5)+2.*py*pz*cov(5,6)+2.*pz*px*cov(4,6))/pr;
  }

  //////////////////////////
  // Main interface function
  //////////////////////////

  void TrkExtDiag::trkExtDiag(const art::Event & evt, const KalRep & krep, const TrkExtTraj & trkext) {
    TrkHitVector const& hots = krep.hitVector();
    _nhots = readHit(evt, hots);
    _ntrks = readTrk(evt, krep);
    _nsim = readMC (evt, krep, hots);
    _next = readExt(evt, trkext);
    _extdiag->Fill();
  }

  void TrkExtDiag::trkExtDiag() {
    _nhots = 0;
    _ntrks = 0;
    _trkl0 = _trkl1 = 0;
    _trk0.d0 = _trk0.z0 = _trk0.p0 = _trk0.om = _trk0.td = 0;
    _trk1.d0 = _trk1.z0 = _trk1.p0 = _trk1.om = _trk1.td = 0;
    _nsim = 0;
    _simp0 = _simt0 = 0;
    _nmcpa = _nmcst = 0;
    _vdsi.status = 0; 
    _vdso.status = 0; 
    _vdtf.status = 0; 
    _vdtm.status = 0; 
    _vdtb.status = 0; 
    _next = 0;
    _nextpa = _nextst = 0;
    _extpadptot =  _extstdptot = 0;
    _extvdsi.status = 0;
    _extvdso.status = 0;
    _extvdtf.status = 0;
    _extvdtm.status = 0;
    _extvdtb.status = 0;
    _extdiag->Fill();
  }



  /////////////
  // readHit //
  /////////////


  unsigned int TrkExtDiag::readHit(const art::Event &evt, const TrkHitVector & hots) {
    unsigned int i = 0;
    for (auto iter = hots.begin() ; iter != hots.end() ; ++iter) {
      const TrkHit * hit  = *iter;
      const mu2e::TrkStrawHit* strawHit = dynamic_cast<const mu2e::TrkStrawHit*>(hit);
      if (strawHit) {
        const HepPoint &point = strawHit->hitTraj()->position(strawHit->hitLen());
        _hotx[i] = point.x();
        _hoty[i] = point.y();
        _hotz[i] = point.z();
        _hott0[i] = strawHit->hitT0()._t0 / CLHEP::ns;
        ++i;
      }
      if (i >= MAXNHOT) {
        cerr << "TrkExtDiag warning : MAXNHOT reached. Remaining are ignored" <<endl;
        return i;
      }
    }
    return i;
  }



  /////////////
  // readTrk //
  /////////////

  unsigned int TrkExtDiag::readTrk(const art::Event &evt, const KalRep & krep){

    unsigned int ntrks = 50-1;
    _trkl0 = krep.startValidRange();
    _trkl1 = krep.endValidRange();
    double step = (_trkl1 - _trkl0)/double(ntrks);

    unsigned int i;
    for( i = 0; i <= ntrks; ++i) {
      double fltlen = _trkl0 + i * step;
      if (i == ntrks) fltlen = _trkl1;
      _trkl[i] = fltlen;
      HepPoint p = krep.position(fltlen);
      _trkx[i] = p.x();
      _trky[i] = p.y();
      _trkz[i] = p.z();
      _trkp[i] = krep.momentum(fltlen).mag();
      _trkpx[i] = krep.momentum(fltlen).x();
      _trkpy[i] = krep.momentum(fltlen).y();
      _trkpz[i] = krep.momentum(fltlen).z();
    }
    HelixParams parl0 = krep.helix(_trkl0);
    _trk0.d0 = parl0.d0();
    _trk0.p0 = parl0.phi0();
    _trk0.om = parl0.omega();
    _trk0.z0 = parl0.z0();
    _trk0.td = parl0.tanDip();
    const HepSymMatrix & cov0 = parl0.covariance();
    _trk0.ed0 = sqrt(cov0(1,1));
    _trk0.ep0 = sqrt(cov0(2,2)); 
    _trk0.eom = sqrt(cov0(3,3));
    _trk0.ez0 = sqrt(cov0(4,4)); 
    _trk0.etd = sqrt(cov0(5,5));
    HelixParams parl1 = krep.helix(_trkl1);
    _trk1.d0 = parl1.d0();
    _trk1.p0 = parl1.phi0();
    _trk1.om = parl1.omega();
    _trk1.z0 = parl1.z0();
    _trk1.td = parl1.tanDip();
    const HepSymMatrix & cov1 = parl1.covariance();
    _trk1.ed0 = sqrt(cov1(1,1));
    _trk1.ep0 = sqrt(cov1(2,2)); 
    _trk1.eom = sqrt(cov1(3,3));
    _trk1.ez0 = sqrt(cov1(4,4)); 
    _trk1.etd = sqrt(cov1(5,5));
    return (unsigned int)i;
  }



  ////////////
  // readMC //
  ////////////

  unsigned int TrkExtDiag::readMC(const art::Event &evt, const KalRep & krep, const TrkHitVector & hot){

    unsigned int i, j;
    bool readflag = false;
    unsigned int minsim = 0;
    _simtp.z = 99999;
    _simtpqual = -1;
    for (auto iter = hot.begin() ; iter != hot.end() ; ++iter) {
      if (readflag) break;
      const TrkHit * hit  = *iter;
      const mu2e::TrkStrawHit* trkstrawHit = dynamic_cast<const mu2e::TrkStrawHit*>(hit);
      if (trkstrawHit == nullptr) {
        cerr << "TrkExtDiag warning : no trkstrawHit" <<endl;
        continue;
      }
      const ComboHit& comboHit = trkstrawHit->comboHit();

      art::Handle<ComboHitCollection> shcHandle;
      evt.getByLabel(_makerModuleLabel, shcHandle);
      if (!(shcHandle.isValid())) {
        cerr << "TrkExtDiag warning : shcHandle invalid" << endl;
        continue;
      }

      ComboHitCollection const& shc = *shcHandle;
      int hitid = -1;

      for (i = 0 ; i <shc.size() ; ++i) {
        ComboHit  const& sh(shc.at(i));
        if (&sh == &comboHit) {
          hitid = int(i);
          break;
        }
      }
      if (hitid <0) {
        cerr << "TrkExtDiag warning : matching straw hit not found" << endl;
        continue;
      }
  
      art::Handle<PtrStepPointMCVectorCollection> stepsCHandle;
      evt.getByLabel(_makerModuleLabel, "StrawHitMCPtr", stepsCHandle);
      if (!(stepsCHandle.isValid())) {
        cerr << "TrkExtDiag warning : stepsCHandle invalid" << endl;
        continue;
      }
      PtrStepPointMCVectorCollection const & stepsC = *stepsCHandle;
  
      if (stepsC.size() <=0) {
        cerr << "TrkExtDiag warning : no stepsC" << endl;
        continue;
      }
  
      PtrStepPointMCVector const & steps (stepsC.at(hitid));
  
      if (steps.size() <=0) {
        cerr << "TrkExtDiag warning : no steps" << endl;
        continue;
      }
  
      StepPointMC const & step = *steps.at(0);
  
      SimParticle const &sim = *step.simParticle();
      
      cet::map_vector_key simid =  sim.id();
  
      art::Handle<PointTrajectoryCollection> trajHandle;
      evt.getByLabel(_g4ModuleLabel, trajHandle);
      if (!(trajHandle.isValid())) {
        cerr << "TrkExtDiag warning : tragHandle invalid" << endl;
        continue;
      }

      const PointTrajectoryCollection & trajC = *trajHandle;
 
      const PointTrajectory* traj = trajC.getOrNull(simid);
      if (traj == nullptr) {
        cerr << "TrkExtDiag warning : traj invalid" << endl;
        continue;
      }
  
      const vector<Hep3Vector>& pvec = traj->points();

      double simstep = 1.;
      if (pvec.size() >= MAXNSIM) {
        simstep = (double)pvec.size() / (double)MAXNSIM;
      }
      minsim = pvec.size();
      if (minsim > MAXNSIM) minsim = MAXNSIM;

      for (i = 0 ; i < minsim ; ++i) {
        j = (int)(simstep*i+0.5);
        if (j >= pvec.size()) j = pvec.size() -1;
        Hep3Vector position = pvec[j];
        _simx[i] = position.x();
        _simy[i] = position.y();
        _simz[i] = position.z() + 1800;
      }
      _simtp.z = 99999;
      _simtpqual = -1;
      
      for (i = 1 ; i < pvec.size()-1 ; ++i) {
        double z1 = pvec[i-1].z();
        double z2 = pvec[i].z();
        double z3 = pvec[i+1].z();
        if (z2+1800.>0) break;
        if ((z3-z2)*(z1-z2)> 0) {
          if (i>10 && i+10<pvec.size()) {
            double z101 = pvec[i-10].z();
            double z102 = pvec[i].z();
            double z103 = pvec[i+10].z();
            if ((z103-z102)*(z101-z102)> 0) _simtpqual = 3;
            else _simtpqual = 2; 
          }
          else {
            _simtpqual = 1;
          }
          double ic = findTurnAroundSim (i, z1, z2, z3);
          if (ic>0) {
            double i1 = (double)((int)ic);
            double i2 = i1+1;
            Hep3Vector p1 = pvec[i1];
            Hep3Vector p2 = pvec[i2];
            _simtp.x = interpolate2(ic, i1, i2,  p1.x(), p2.x());
            _simtp.y = interpolate2(ic, i1, i2,  p1.y(), p2.y());
            _simtp.z = interpolate2(ic, i1, i2,  p1.z(), p2.z()) + 1800.;
          }
          break;
        }
      }
  
      _simt0 = sim.startGlobalTime();
      _simp0 = sim.startMomentum().vect().mag();

      readflag = true;

      //readvd
      _vdsi.clear();
      _vdso.clear();
      _vdtf.clear();
      _vdtm.clear();
      _vdtb.clear();
      TrkExtMCHits mcvd(evt, _g4ModuleLabel, "virtualdetector", simid);
      int hepid = sim.pdgId();
      int sign;
      if (hepid == 11 || hepid == 13 || hepid == -211 || hepid == -321 || hepid == -2212) sign = 1;
      else if (hepid == -11 || hepid == -13 || hepid == 211 || hepid == 321 || hepid == 2212) sign = -1;
      else sign = 0;
      for (i = 0 ; i <mcvd.getNClusters() ; ++i) {
        const StepPointMC & firsthit = (mcvd.getCluster(i))[0];
        Hep3Vector mom = mcvd.momentum(i);
        Hep3Vector pos = mcvd.position(i) - _origin;
        HepVector hpar;
        unsigned int index;
        if (mom.z() >=0) index = 0;
        else index = 1;
        switch (firsthit.volumeId()) {
          case VirtualDetectorId::ST_In :
            if ( _vdsi.status & (unsigned int)(index+1) ) {
              cerr << "TrkExtDiag warning : Too many VD hit at " << firsthit.volumeId() << ". Ignored." << endl;
              continue;
            }
            _vdsi.x[index] = pos.x();
            _vdsi.y[index] = pos.y();
            _vdsi.z[index] = pos.z();
            _vdsi.px[index] = mom.x();
            _vdsi.py[index] = mom.y();
            _vdsi.pz[index] = mom.z();
            _vdsi.p[index] = mom.mag();
            hpar = getHelixParameters (pos, mom, sign);
            _vdsi.d0[index] = hpar(1);
            _vdsi.p0[index] = hpar(2);
            _vdsi.om[index] = hpar(3);
            _vdsi.z0[index] = hpar(4);
            _vdsi.td[index] = hpar(5);
            _vdsi.status = _vdsi.status | (unsigned int)(index+1);
            break;
          case VirtualDetectorId::ST_Out :
            if ( _vdso.status & (unsigned int)(index+1) ) {
              cerr << "TrkExtDiag warning : Too many VD hit at " << firsthit.volumeId() << ". Ignored." << endl;
              continue;
            }
            _vdso.x[index] = pos.x();
            _vdso.y[index] = pos.y();
            _vdso.z[index] = pos.z();
            _vdso.px[index] = mom.x();
            _vdso.py[index] = mom.y();
            _vdso.pz[index] = mom.z();
            _vdso.p[index] = mom.mag();
            hpar = getHelixParameters (pos, mom, sign);
            _vdso.d0[index] = hpar(1);
            _vdso.p0[index] = hpar(2);
            _vdso.om[index] = hpar(3);
            _vdso.z0[index] = hpar(4);
            _vdso.td[index] = hpar(5);
            _vdso.status = _vdso.status | (unsigned int)(index+1);
            break;
          case VirtualDetectorId::TT_FrontHollow :
          case VirtualDetectorId::TT_FrontPA :
            if ( _vdtf.status & (unsigned int)(index+1) ) {
              cerr << "TrkExtDiag warning : Too many VD hit at " << firsthit.volumeId() << ". Ignored." << endl;
              continue;
            }
            _vdtf.x[index] = pos.x();
            _vdtf.y[index] = pos.y();
            _vdtf.z[index] = pos.z();
            _vdtf.px[index] = mom.x();
            _vdtf.py[index] = mom.y();
            _vdtf.pz[index] = mom.z();
            _vdtf.p[index] = mom.mag();
            hpar = getHelixParameters (pos, mom, sign);
            _vdtf.d0[index] = hpar(1);
            _vdtf.p0[index] = hpar(2);
            _vdtf.om[index] = hpar(3);
            _vdtf.z0[index] = hpar(4);
            _vdtf.td[index] = hpar(5);
            _vdtf.status = _vdtf.status | (unsigned int)(index+1);
            break;
          case VirtualDetectorId::TT_MidInner :
          case VirtualDetectorId::TT_Mid :
            if ( _vdtm.status & (unsigned int)(index+1) ) {
              cerr << "TrkExtDiag warning : Too many VD hit at " << firsthit.volumeId() << ". Ignored." << endl;
              continue;
            }
            _vdtm.x[index] = pos.x();
            _vdtm.y[index] = pos.y();
            _vdtm.z[index] = pos.z();
            _vdtm.px[index] = mom.x();
            _vdtm.py[index] = mom.y();
            _vdtm.pz[index] = mom.z();
            _vdtm.p[index] = mom.mag();
            hpar = getHelixParameters (pos, mom, sign);
            _vdtm.d0[index] = hpar(1);
            _vdtm.p0[index] = hpar(2);
            _vdtm.om[index] = hpar(3);
            _vdtm.z0[index] = hpar(4);
            _vdtm.td[index] = hpar(5);
            _vdtm.status = _vdtm.status | (unsigned int)(index+1);
            break;
          case VirtualDetectorId::TT_Back :
            if ( _vdtb.status & (unsigned int)(index+1) ) {
              cerr << "TrkExtDiag warning : Too many VD hit at " << firsthit.volumeId() << ". Ignored." << endl;
              continue;
            }
            _vdtb.x[index] = pos.x();
            _vdtb.y[index] = pos.y();
            _vdtb.z[index] = pos.z();
            _vdtb.px[index] = mom.x();
            _vdtb.py[index] = mom.y();
            _vdtb.pz[index] = mom.z();
            _vdtb.p[index] = mom.mag();
            hpar = getHelixParameters (pos, mom, sign);
            _vdtb.d0[index] = hpar(1);
            _vdtb.p0[index] = hpar(2);
            _vdtb.om[index] = hpar(3);
            _vdtb.z0[index] = hpar(4);
            _vdtb.td[index] = hpar(5);
            _vdtb.status = _vdtb.status | (unsigned int)(index+1);
            break;
          default :
            break;
        }
      }

      // pa 
      TrkExtMCHits mcpa(evt, _g4ModuleLabel, "protonabsorber", simid);
      _nmcpa = mcpa.getNClusters();
      for (i = 0 ; i < _nmcpa ; ++i) {
        if (i >= MAXNPA) {
          cerr << "ReadTrkExt : too many PA hits " << mcpa.getNClusters() << endl;
          break;
        }
        Hep3Vector mom = mcpa.momentum(i);
        Hep3Vector pos = mcpa.position(i) - _origin;
        _mcpapx[i] = mom.x();
        _mcpapy[i] = mom.y();
        _mcpapz[i] = mom.z();
        _mcpap[i] = mom.mag();
        _mcpaz[i] = pos.z();
        _mcpadp[i] = mcpa.deltap(i);
        _mcpade[i] = mcpa.eDep(i);
        _mcpadei[i] = mcpa.ionizingEdep(i);
        _mcpadeni[i] = mcpa.nonIonizingEdep(i);
      }

      // st
      TrkExtMCHits mcst(evt, _g4ModuleLabel, "stoppingtarget", simid);
      _nmcst = mcst.getNClusters();
      for (i = 0 ; i < _nmcst ; ++i) {
        if (i >= MAXNST) {
          cerr << "ReadTrkExt : too many ST hits " << mcst.getNClusters() << endl;
          break;
        }
        Hep3Vector mom = mcst.momentum(i);
        Hep3Vector pos = mcst.position(i) - _origin;
        _mcstpx[i] = mom.x();
        _mcstpy[i] = mom.y();
        _mcstpz[i] = mom.z();
        _mcstp[i] = mom.mag();
        _mcstx[i] = pos.x();
        _mcsty[i] = pos.y();
        _mcstz[i] = pos.z();
        _mcstt[i] = mcst.time(i);
        _mcstdp[i] = mcst.deltap(i);
        _mcstde[i] = mcst.eDep(i);
        _mcstdei[i] = mcst.ionizingEdep(i);
        _mcstdeni[i] = mcst.nonIonizingEdep(i);
      }

    } // end of hot loop
  
    if (!readflag) {
      cerr << "TrkExtDiag Error: Cannot read MC information" << endl;
      return 0;
    }
    return minsim;
  }



  /////////////
  // readExt //
  /////////////

  unsigned int TrkExtDiag::readExt(const art::Event &evt, const TrkExtTraj & trkext){
    unsigned int i;
    // trkext 
    double x, y, z, px, py, pz, p, rho;
    int turnaround = -1;
    for (i = 0 ; i < trkext.size() ; ++i) {
      const TrkExtTrajPoint & hit = trkext[i];
      _extx[i] = x = hit.x();
      _exty[i] = y = hit.y();
      _extz[i] = z = hit.z();
      _extrho[i] = rho = hit.rho();
      _extp[i] = p = hit.p();
      _extpx[i] = px = hit.px();
      _extpy[i] = py = hit.py();
      _extpz[i] = pz = hit.pz();
      _extex[i] = hit.ex(); 
      _extey[i] = hit.ey(); 
      _extez[i] = hit.ez(); 
      _extepx[i] = hit.epx();
      _extepy[i] = hit.epy();
      _extepz[i] = hit.epz();
      _extep[i] = safeSqrt( hit.covpxpx()*px*px + hit.covpypy()*py*py + hit.covpzpz()*pz*pz + 2.*px*py*hit.covpxpy() + 2.*py*pz*hit.covpypz() + 2.*pz*px*hit.covpxpz() ) / p;
      _exter[i] = safeSqrt( hit.covxx()*x*x + hit.covyy()*y*y + 2.*x*y*hit.covxy() ) / rho;
      _extvid[i] = hit.volumeId();
      _extt[i] = hit.flightTime();
      _exts[i] = hit.flightLength();
      if (i>0) {
        if (_extpz[i] * _extpz[i-1] <=0) turnaround = i;
      }
    }

    // turn around point
    _extitp = turnaround;
    _exttp.s = 0;
    if (int(_extitp-1) >=0 && int(_extitp +1)<int(trkext.size())) {
      int i1 = _extitp-1;
      int i2 = _extitp;
      int i3 = _extitp+1;
      double s1 = _exts[i1];
      double s2 = _exts[i2];
      double s3 = _exts[i3];
      double s = findTurnAround (s1, s2, s3, _extpz[i1], _extpz[i2], _extpz[i3]);
      if (s !=0) {
        _exttp.x =  interpolate(s, s1, s2, s3, _extx[i1], _extx[i2], _extx[i3]);
        _exttp.y =  interpolate(s, s1, s2, s3, _exty[i1], _exty[i2], _exty[i3]);
        _exttp.z =  interpolate(s, s1, s2, s3, _extz[i1], _extz[i2], _extz[i3]);
        _exttp.px = interpolate(s, s1, s2, s3, _extpx[i1], _extpx[i2], _extpx[i3]);
        _exttp.py = interpolate(s, s1, s2, s3, _extpy[i1], _extpy[i2], _extpy[i3]);
        _exttp.pz = 0;
        _exttp.s =  s;
        _exttp.t =  interpolate(s, s1, s2, s3, _extt[i1], _extt[i2], _extt[i3]);
      }
    }

    // pa and st
    _nextpa = trkext.getNPAHits();
    for ( i = 0 ; i < _nextpa ; ++i) {
      _extpaz[i] = trkext.getMeanPAHit(i).z();
      _extpadp[i] = trkext.getDeltapPA(i);
    }
    _nextst = trkext.getNSTHits();
    for ( i = 0 ; i < _nextst ; ++i) {
      TrkExtTrajPoint meansthit = trkext.getMeanSTHit(i);
      _extstx[i] = meansthit.x();
      _extsty[i] = meansthit.y();
      _extstz[i] = meansthit.z();
      _extstt[i] = meansthit.flightTime();
      _extstdp[i] = trkext.getDeltapST(i);
    }
    _extpadptot = trkext.getDeltapPA();
    _extstdptot = trkext.getDeltapST();

    // vd
    _extvdsi.clear();
    _extvdso.clear();
    _extvdtf.clear();
    _extvdtm.clear();
    _extvdtb.clear();
    int hepid = trkext.hepid();
    int sign;
    if (hepid == 11 || hepid == 13 || hepid == -211 || hepid == -321 || hepid == -2212) sign = 1;
    else if (hepid == -11 || hepid == -13 || hepid == 211 || hepid == 321 || hepid == 2212) sign = -1;
    else sign = 0;

    std::vector<TrkExtTrajPoint> vdsihits = trkext.getPointsAtZ(_vdzsi);
    for (i = 0 ; i < vdsihits.size() ; ++i) {
      unsigned int index;
      if (vdsihits[i].pz()>=0) index = 0;
      else index = 1;
      if (_extvdsi.status & (unsigned int)(index+1)) {
        cout << "TrkExtDiag warning : more than 2 same directional hits in VDsi : " << endl;
        continue;
      }
      _extvdsi.r[index] = vdsihits[i].rho();
      _extvdsi.p[index] = vdsihits[i].p();
      const Hep3Vector & pos = vdsihits[i].position();
      const Hep3Vector & mom = vdsihits[i].momentum();
      const HepMatrix & cov = vdsihits[i].covariance();
      _extvdsi.er[index] = getRadialError(pos, cov);
      _extvdsi.ep[index] = getMomentumError(mom, cov);
      //cout << _extvdsi.er[index] <<", " <<_extvdsi.ep[index] << endl;
      HepVector hpar = getHelixParameters (pos, mom, sign);
      Hep3Vector poserr = vdsihits[i].positionError();
      Hep3Vector momerr = vdsihits[i].momentumError();
      HepVector hepar = getHelixParametersErr (pos, mom, poserr, momerr, sign);
      _extvdsi.d0[index] = hpar(1);
      _extvdsi.p0[index] = hpar(2);
      _extvdsi.om[index] = hpar(3);
      _extvdsi.z0[index] = hpar(4);
      _extvdsi.td[index] = hpar(5);
      _extvdsi.ed0[index] = hepar(1);
      _extvdsi.ep0[index] = hepar(2);
      _extvdsi.eom[index] = hepar(3);
      _extvdsi.ez0[index] = hepar(4);
      _extvdsi.etd[index] = hepar(5);
      _extvdsi.status = _extvdsi.status | (unsigned int)(index+1);
    }

    std::vector<TrkExtTrajPoint> vdsohits = trkext.getPointsAtZ(_vdzso);
    for (i = 0 ; i < vdsohits.size() ; ++i) {
      unsigned int index;
      if (vdsohits[i].pz()>=0) index = 0;
      else index = 1;
      if (_extvdso.status & (unsigned int)(index+1)) {
        cout << "TrkExtDiag warning : more than 2 same directional hits in VDsi : " << endl;
        continue;
      }
      _extvdso.r[index] = vdsohits[i].rho();
      _extvdso.p[index] = vdsohits[i].p();
      const Hep3Vector & pos = vdsohits[i].position();
      const Hep3Vector & mom = vdsohits[i].momentum();
      const HepMatrix & cov = vdsohits[i].covariance();
      _extvdso.er[index] = getRadialError(pos, cov);
      _extvdso.ep[index] = getMomentumError(mom, cov);
      HepVector hpar = getHelixParameters (pos, mom, sign);
      Hep3Vector poserr = vdsohits[i].positionError();
      Hep3Vector momerr = vdsohits[i].momentumError();
      HepVector hepar = getHelixParametersErr (pos, mom, poserr, momerr, sign);
      _extvdso.d0[index] = hpar(1);
      _extvdso.p0[index] = hpar(2);
      _extvdso.om[index] = hpar(3);
      _extvdso.z0[index] = hpar(4);
      _extvdso.td[index] = hpar(5);
      _extvdso.ed0[index] = hepar(1);
      _extvdso.ep0[index] = hepar(2);
      _extvdso.eom[index] = hepar(3);
      _extvdso.ez0[index] = hepar(4);
      _extvdso.etd[index] = hepar(5);
      _extvdso.status = _extvdso.status | (unsigned int)(index+1);
    }

    std::vector<TrkExtTrajPoint> vdtfhits = trkext.getPointsAtZ(_vdztf);
    for (i = 0 ; i < vdtfhits.size() ; ++i) {
      unsigned int index;
      if (vdtfhits[i].pz()>=0) index = 0;
      else index = 1;
      if (_extvdtf.status & (unsigned int)(index+1)) {
        cout << "TrkExtDiag warning : more than 2 same directional hits in VDsi : " << endl;
        continue;
      }
      _extvdtf.r[index] = vdtfhits[i].rho();
      _extvdtf.p[index] = vdtfhits[i].p();
      const Hep3Vector & pos = vdtfhits[i].position();
      const Hep3Vector & mom = vdtfhits[i].momentum();
      const HepMatrix & cov = vdtfhits[i].covariance();
      _extvdtf.er[index] = getRadialError(pos, cov);
      _extvdtf.ep[index] = getMomentumError(mom, cov);
      HepVector hpar = getHelixParameters (pos, mom, sign);
      Hep3Vector poserr = vdtfhits[i].positionError();
      Hep3Vector momerr = vdtfhits[i].momentumError();
      HepVector hepar = getHelixParametersErr (pos, mom, poserr, momerr, sign);
      _extvdtf.d0[index] = hpar(1);
      _extvdtf.p0[index] = hpar(2);
      _extvdtf.om[index] = hpar(3);
      _extvdtf.z0[index] = hpar(4);
      _extvdtf.td[index] = hpar(5);
      _extvdtf.ed0[index] = hepar(1);
      _extvdtf.ep0[index] = hepar(2);
      _extvdtf.eom[index] = hepar(3);
      _extvdtf.ez0[index] = hepar(4);
      _extvdtf.etd[index] = hepar(5);
      _extvdtf.status = _extvdtf.status | (unsigned int)(index+1);
    }

    std::vector<TrkExtTrajPoint> vdtmhits = trkext.getPointsAtZ(_vdztm);
    for (i = 0 ; i < vdtmhits.size() ; ++i) {
      unsigned int index;
      if (vdtmhits[i].pz()>=0) index = 0;
      else index = 1;
      if (_extvdtm.status & (unsigned int)(index+1)) {
        cout << "TrkExtDiag warning : more than 2 same directional hits in VDsi : " << endl;
        continue;
      }
      _extvdtm.r[index] = vdtmhits[i].rho();
      _extvdtm.p[index] = vdtmhits[i].p();
      const Hep3Vector & pos = vdtmhits[i].position();
      const Hep3Vector & mom = vdtmhits[i].momentum();
      const HepMatrix & cov = vdtmhits[i].covariance();
      _extvdtm.er[index] = getRadialError(pos, cov);
      _extvdtm.ep[index] = getMomentumError(mom, cov);
      HepVector hpar = getHelixParameters (pos, mom, sign);
      Hep3Vector poserr = vdtmhits[i].positionError();
      Hep3Vector momerr = vdtmhits[i].momentumError();
      HepVector hepar = getHelixParametersErr (pos, mom, poserr, momerr, sign);
      _extvdtm.d0[index] = hpar(1);
      _extvdtm.p0[index] = hpar(2);
      _extvdtm.om[index] = hpar(3);
      _extvdtm.z0[index] = hpar(4);
      _extvdtm.td[index] = hpar(5);
      _extvdtm.ed0[index] = hepar(1);
      _extvdtm.ep0[index] = hepar(2);
      _extvdtm.eom[index] = hepar(3);
      _extvdtm.ez0[index] = hepar(4);
      _extvdtm.etd[index] = hepar(5);
      _extvdtm.status = _extvdtm.status | (unsigned int)(index+1);
    }

    std::vector<TrkExtTrajPoint> vdtbhits = trkext.getPointsAtZ(_vdztb);
    for (i = 0 ; i < vdtbhits.size() ; ++i) {
      unsigned int index;
      if (vdtbhits[i].pz()>=0) index = 0;
      else index = 1;
      if (_extvdtb.status & (unsigned int)(index+1)) {
        cout << "TrkExtDiag warning : more than 2 same directional hits in VDsi : " << endl;
        continue;
      }
      _extvdtb.r[index] = vdtbhits[i].rho();
      _extvdtb.p[index] = vdtbhits[i].p();
      const Hep3Vector & pos = vdtbhits[i].position();
      const Hep3Vector & mom = vdtbhits[i].momentum();
      const HepMatrix & cov = vdtbhits[i].covariance();
      _extvdtb.er[index] = getRadialError(pos, cov);
      _extvdtb.ep[index] = getMomentumError(mom, cov);
      HepVector hpar = getHelixParameters (pos, mom, sign);
      Hep3Vector poserr = vdtbhits[i].positionError();
      Hep3Vector momerr = vdtbhits[i].momentumError();
      HepVector hepar = getHelixParametersErr (pos, mom, poserr, momerr, sign);
      _extvdtb.d0[index] = hpar(1);
      _extvdtb.p0[index] = hpar(2);
      _extvdtb.om[index] = hpar(3);
      _extvdtb.z0[index] = hpar(4);
      _extvdtb.td[index] = hpar(5);
      _extvdtb.ed0[index] = hepar(1);
      _extvdtb.ep0[index] = hepar(2);
      _extvdtb.eom[index] = hepar(3);
      _extvdtb.ez0[index] = hepar(4);
      _extvdtb.etd[index] = hepar(5);
      _extvdtb.status = _extvdtb.status | (unsigned int)(index+1);
    }

    // exitcode
    _exitcode = trkext.exitCode();

    return trkext.size();
  }

}
