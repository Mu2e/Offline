// Generate anti-proton with a shape
//
// Ying Wang, Zhengyun You, 2018-09-13

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"

#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"


namespace mu2e {

  //================================================================
  namespace {
    struct VDHit {
      float x;
      float y;
      float z;
      float time;
      float px;
      float py;
      float pz;
      int   pdgId;

      VDHit() : x(std::numeric_limits<double>::quiet_NaN())
              , y(std::numeric_limits<double>::quiet_NaN())
              , z(std::numeric_limits<double>::quiet_NaN())
              , time(std::numeric_limits<double>::quiet_NaN())
              , px(std::numeric_limits<double>::quiet_NaN())
              , py(std::numeric_limits<double>::quiet_NaN())
              , pz(std::numeric_limits<double>::quiet_NaN())
              , pdgId(0)
      {}

    }; // struct VDHit
  } // namespace {}

  //================================================================
  class AntiProtonGun : public art::EDProducer {

    PDGCode::type pdgId_;
    int verbosityLevel_;

    static double zLambda(double *val, double *par);
    static double plmax(double theta, double xx);
    static double dsigma(double *val, double *parm);

    bool rotateTarget_;

    double phiMin_;   // rotate an angle between phiMin_ and phiMax_
    double phiMax_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randFlat_;
    CLHEP::RandGaussQ randGaussQ_;

    CLHEP::HepRotation targetRotation_; // rotates target frame to Mu2e frame
    CLHEP::Hep3Vector  targetCenter_;
  
    TF1 *f1;
    TF2 *fsig;

    bool firstEvent_;

  public:
    explicit AntiProtonGun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
    virtual void beginRun(art::Run& run);
  };

  //================================================================
  AntiProtonGun::AntiProtonGun(const fhicl::ParameterSet& pset)
    : art::EDProducer{pset}
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , rotateTarget_(pset.get<bool>("rotateTarget", false))
    , phiMin_(pset.get<double>("phiMin", 0.))
    , phiMax_(pset.get<double>("phiMax", 360.))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randFlat_(eng_)
    , randGaussQ_(eng_, 0, 1.0)
    , firstEvent_(true)
  {
    produces<mu2e::GenParticleCollection>();

    if(verbosityLevel_ > 0) {
      std::cout<<"AntiProtonGun: producing particle "
               <<"..."
               <<std::endl;
    }

    phiMin_   *= CLHEP::deg;
    phiMax_   *= CLHEP::deg;
  }

  //================================================================
  void AntiProtonGun::beginRun(art::Run& run) {

    if(verbosityLevel_ > 0) {
      std::cout<<"AntiProtonGun beginRun "
               <<std::endl;
    }
  }

  //================================================================
  double AntiProtonGun::zLambda(double *val, double *par) {

    par = 0;
    double kLambda = 99.4;
    double kTgtHL = 80.0;
    double value=1/kLambda*exp(-(kTgtHL-val[0])/kLambda);

    return value;
  }


  //================================================================
  double AntiProtonGun::plmax(double theta, double xx) {

    double pLab=8.8;
    double x=xx;
    double Mp=0.938272029;
    double plmax0 = 0;
    double smx = 3.*Mp;
    double e0 = sqrt(pLab*pLab + Mp*Mp);  // beam energy
    double s = 2*x*Mp*e0 + x*x*Mp*Mp + Mp*Mp;
    double rootS = sqrt(s);
    double eMax = (s + Mp*Mp - smx*smx)/(2.*rootS);

    double gacm = (e0 + x*Mp)/rootS;  // = rootS/2/mP =2.28387 pA
    double gabcm = pLab/(e0+x*Mp);  // pA

    double a = gacm*gacm*gabcm*gabcm*cos(theta)*cos(theta)-gacm*gacm;
    double b = 2*eMax*gabcm*gacm*cos(theta);
    double c =eMax*eMax-gacm*gacm*Mp*Mp;
    double disc = b*b-4*a*c;

    if (disc >= 0) {
        plmax0 = (-b-sqrt(disc))/(2*a);
    }

    return plmax0;
  } // plmax

  //================================================================
  double AntiProtonGun::dsigma(double *val, double *parm) {

    parm=0;
    double x=1.9;
    double total=0;
    double theta=val[0];
    double pbarmax=plmax(theta, x);
    double pbar=val[1];

    double pLab=8.8;
    double Mp=0.938272029;
 
    double kpar[10]={
    0.1699,  //___0
    10.28,   //___1
    2.269,   //___2
    3.707,   //___3
    0.009,   //___4
    0.4812,  //___5
    3.3600,  //___6
    0.06394, //___7
    -0.1824, //___8
    2.4850   //___9
    };

    double smx = 3.*Mp;
    double e0 = sqrt(pLab*pLab + Mp*Mp);  // beam energy
    double s = 2*x*Mp*e0 + x*x*Mp*Mp + Mp*Mp;
    double rootS = sqrt(s);
    double eMax = (s + Mp*Mp - smx*smx)/(2.*rootS);
    double gacm = (e0 + x*Mp)/rootS;  // = rootS/2/mP =2.28387 pA
    double gabcm = pLab/(e0+x*Mp);  // pA
    double E_cm=gacm*(sqrt(pbar*pbar+Mp*Mp)-gabcm*pbar*cos(theta));

    if (E_cm < eMax && pbar<=pbarmax) {
      double x_R=E_cm/eMax;

      double pt=pbar*sin(theta);
      double Einc=(s/2./Mp-2*Mp)*1000;

      double T1=pow(184,kpar[0]*log(sqrt(s)/kpar[1])*pt);
      double T2=pow(1-x_R,kpar[2]*log(sqrt(s)));
      double T_R=exp(-kpar[3]*x_R);
      double T4=kpar[4]*pow(sqrt(s),kpar[5])*exp(-kpar[6]*pt);
      double T5=kpar[7]*pow(sqrt(s),kpar[8])*exp(-kpar[9]*pt);
      double factor=50;
      double sig0=45*pow(184,0.7)*(1+0.016*sin(5.3-2.63*log(184))); //wt 184
      double sig_in=factor*sig0*(1-0.62*exp(-Einc/200.)*sin(10.9*pow(Einc,-0.28)));

      total=sig_in*T1*T2*T_R*(T4+T5);
    }

    return total;
  } // dsigma

  //================================================================
  void AntiProtonGun::produce(art::Event& event) {

    if (firstEvent_) {
      GeomHandle<ProductionTarget> target;
      targetRotation_ = target->protonBeamRotation();
      targetCenter_ = target->position();
      f1 = new TF1("f1", AntiProtonGun::zLambda, -80.0, 80.0);
      fsig = new TF2("fsig", AntiProtonGun::dsigma, 0, M_PI, 0, 10.0);  // kPbarLMax=10.0;
      fsig->SetNpx(180);   // number of cells the x axis divided into while randomizing
      fsig->SetNpy(1000);  // number of cells the y axis divided into

      if(verbosityLevel_ > 1) {
        std::cout << "targetCenter@  " << targetCenter_ << std::endl;
        std::cout << "targetRotation " << targetRotation_ << std::endl;
      }
      firstEvent_ = false;
    }

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    //const CLHEP::Hep3Vector pos(3904, 0, -6164.5);
    //const CLHEP::Hep3Vector mom(-241.922, 0, -970.296);

    double phi = 2. * M_PI * randFlat_.fire();
    if(verbosityLevel_ > 2) std::cout << "phi " << phi << std::endl;
    double theta;
    double pbarL;
    fsig->GetRandom2(theta, pbarL);
    CLHEP::Hep3Vector mom_pbarTarget(-pbarL*sin(theta)*cos(phi), -pbarL*sin(theta)*sin(phi), -pbarL*cos(theta));

    double r=1e9;
    while (r>=3.15) { //kTgtR=3.15mm
      r = fabs(randGaussQ_.fire());
    }
    phi = 2. * M_PI * randFlat_.fire(); // position phi does not have to be momentum phi
    double z=f1->GetRandom();
    const double time=(80.0-z)/300.0; //kLS=300.0mm/ns;
    CLHEP::Hep3Vector pos_pbarTarget(r*cos(phi), r*sin(phi), z);

    CLHEP::Hep3Vector mom = targetRotation_*mom_pbarTarget*1000.0;
    CLHEP::Hep3Vector pos = targetRotation_*pos_pbarTarget+targetCenter_;

    CLHEP::Hep3Vector mom_new = mom;
    CLHEP::Hep3Vector pos_new = pos;

    if (rotateTarget_) {
      if(verbosityLevel_ > 1) std::cout << "pos " << pos << " mom " << mom << std::endl;

      CLHEP::Hep3Vector mom_tgt = targetRotation_.inverse()*mom;
      CLHEP::Hep3Vector pos_tgt = targetRotation_.inverse()*(pos-targetCenter_);
      if(verbosityLevel_ > 1) std::cout << "pos_tgt " << pos_tgt << " mom_tgt " << mom_tgt << std::endl;

      // rotate by an random angle around phi and get new p direction in target frame
      double rotatePhi = randFlat_.fire()*(phiMax_-phiMin_);
      if(verbosityLevel_ > 1) std::cout << "rotatePhi " << rotatePhi << std::endl;

      mom_tgt.rotateZ(rotatePhi);
      pos_tgt.rotateZ(rotatePhi);
      if(verbosityLevel_ > 1) std::cout << "pos_tgt_rot " << pos_tgt << " mom_tgt_rot " << mom_tgt << std::endl;

      // rotate around y back to mu2e frame
      mom_new = targetRotation_*mom_tgt;
      pos_new = targetRotation_*pos_tgt+targetCenter_;
      if(verbosityLevel_ > 1) std::cout << "pos_rot " << pos_new << " mom_rot " << mom_new << std::endl;
    }

    if(verbosityLevel_ > 2) std::cout << "pos " << pos_new << " mom " << mom_new << std::endl;

    PDGCode::type pdgId = static_cast<PDGCode::type>(-2212);
    const double mass = GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId).ref().mass().value();
    const double energy = sqrt(mom_new.mag2() + mass*mass);
    CLHEP::HepLorentzVector fourmom(mom_new, energy);

    GenParticle outGen(pdgId, GenId::fromSimParticleCompact, pos_new, fourmom, time);
    output->push_back(outGen);

    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::AntiProtonGun);
