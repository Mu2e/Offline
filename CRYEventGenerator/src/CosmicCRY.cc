// C++ includes.
#include <cmath>
#include <iostream>
#include <cstdlib>

// Framework includes.
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/rm48.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

#include "CRYGenerator.h"
#include "CRYSetup.h"
#include "CRYEventGenerator/inc/CosmicCRY.hh"
// ROOT includes
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

// From CLHEP
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
using CLHEP::RandFlat;
using CLHEP::GeV;

namespace mu2e 
{
  // The following part is needed for the RNG
  template<class T> class RNGWrapper {
    public:
      static void set(T* object, double (T::*func)(void));
      static double rng(void);
    private:
      static T* m_obj;
      static double (T::*m_func)(void);
  };// end of RNGWrapper class

  template<class T> T* RNGWrapper<T>::m_obj;
  template<class T> double (T::*RNGWrapper<T>::m_func)(void);
  template<class T> void RNGWrapper<T>::set(T* object, double (T::*func)(void)) {
    m_obj = object; m_func = func;
  }
  template<class T> double RNGWrapper<T>::rng(void) { return (m_obj->*m_func)(); }


  CosmicCRY::CosmicCRY( art::Run& run,
      const SimpleConfig& config, CLHEP::HepRandomEngine & engine )
    : _verbose(config.getInt("cosmicCRY.verbose", 0) )
      , _doHistograms(config.getBool("cosmicCRY.doHistograms", true) )
      , _saveTree(config.getBool("cosmicCRY.saveTree", false) )
      , _hXZ(NULL)
      , _hY(NULL)
      // , _hPlane(NULL)
      , _hE(NULL)
      , _hTheta(NULL)
      , _hPhi(NULL)
      , _hPtot(NULL)
      , _hPyOverPtot(NULL)
      , _hTime(NULL)
      , _hLiveTime(NULL)
      , _hNegMuKE(NULL)
      , _hPosMuKE(NULL)
      , _hPtypeKE(NULL)
      , _hNSecondaries(NULL)
      , _hSecondPtotVsPrimKE(NULL)
      , _hShowerRadiusVsPrimKE(NULL)
      , _hNSecondariesVsPrimKE(NULL)

      , _muEMin    ( config.getDouble("cosmicCRY.muEMin") )   //in MeV
        , _muEMax    ( config.getDouble("cosmicCRY.muEMax") )   //in MeV
        , _muCosThMin( config.getDouble("cosmicCRY.muCosThMin") )
        , _muCosThMax( config.getDouble("cosmicCRY.muCosThMax") )
        , _muPhiMin(0)
        , _muPhiMax(2.0*M_PI)

        , _returnMuons( config.getBool("cosmicCRY.returnMuons", true))
        , _returnNeutrons( config.getBool("cosmicCRY.returnNeutrons", true))
        , _returnProtons( config.getBool("cosmicCRY.returnProtons", true))
        , _returnGammas( config.getBool("cosmicCRY.returnGammas", true))
        , _returnElectrons( config.getBool("cosmicCRY.returnElectrons", true))
        , _returnPions( config.getBool("cosmicCRY.returnPions", true))
        , _returnKaons( config.getBool("cosmicCRY.returnKaons", true))

        , _month(config.getInt("cosmicCRY.month", 6))
        , _day(config.getInt("cosmicCRY.day", 21))
        , _year(config.getInt("cosmicCRY.year", 2020))
        , _latitude( config.getDouble("cosmicCRY.latitude", 41.8))
        , _altitude( config.getInt("cosmicCRY.altitude", 0))
        , _subboxLength( config.getDouble("cosmicCRY.subboxLength", 100.))

        , _setupString("")
        , _refY0(config.getDouble("cosmicCRY.refY0", 10000.))
        , _refPointChoice(config.getString("cosmicCRY.refPoint", "UNDEFINED"))
        , _directionChoice(config.getString("cosmicCRY.directionChoice", "ALL"))
        , _cosmicReferencePointInMu2e()
        , _vertical(false)
        , _projectToEnvelope(config.getBool("cosmicCRY.projectToEnvelope", false))
        , _geomInfoObtained(false)

        {
          mf::LogInfo log("CosmicCRY");

          _cryDataPath = std::string(std::getenv("CRYDATAPATH"));
          if( _cryDataPath.length() == 0) 
          {
            mf::LogError("CosmicCRY") << "no variable CRYDATAPATH set. Exit now.";
            exit(0);
          }

          createSetupString();
          _crySetup = new CRYSetup(_setupString, _cryDataPath);

          RNGWrapper<CLHEP::HepRandomEngine>::set(
              &engine,
              &CLHEP::HepRandomEngine::flat);
          _crySetup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);

          if (_verbose > 1) 
            engine.showStatus();

          _cryGen = std::make_shared<CRYGenerator>(_crySetup);

          // Histograming
          if (_doHistograms) {
            art::ServiceHandle<art::TFileService> tfs;
            art::TFileDirectory tfdir = tfs->mkdir("CosmicCRY");
            _hXZ = tfdir.make<TH2D>("XZ", "XZ", 
                400, -2.0e5,  2.0e5, 400, -2.0e5, 2.0e5 );
            _hY = tfdir.make<TH1D>("Y", "Y", 100, -1.0e3, 1.0e3 );

            // _hPlane = tfdir.make<TH1D>("Plane", "Plane", 5, 0, 5);
            _hE = tfdir.make<TH1D>("E", "E", 4000, 0, _muEMax);
            _hTheta = tfdir.make<TH1D>("Theta", "Theta",
                200, -M_PI - 0.5, M_PI + 0.5);
            _hPhi = tfdir.make<TH1D>("Phi", "Phi",
                200, -M_PI - 0.5, M_PI + 0.5);

            _hPtot = tfdir.make<TH1D>("Ptot", "Total momentum", 500, 0., 2E6);
            _hPyOverPtot = tfdir.make<TH1D>("PyOverPtot", "Py/Ptot", 200, -20., 20.);
            _hTime = tfdir.make<TH1D>("Time", "Timing of secondary particles",
                1000, -1.E-3, 1E0);
            _hLiveTime = tfdir.make<TH1D>("Livetime", "Live time", 1000, -1.E-3, 1E0);
            _hNegMuKE = tfdir.make<TH1D>("NegMuKE", "Momenta of #mu^-", 4000, 0., 2E6);
            _hPosMuKE = tfdir.make<TH1D>("PosMuKE", "Momenta of #mu^+", 4000, 0., 2E6);
            _hNSecondaries = tfdir.make<TH1D>("nSecondaries", "Number of secondaries",
                1000, 0, 1000);

            _hPtypeKE = tfdir.make<TH2D>("PtypeKE", "Particle type vs energy",
                2000, 0., 2E6, 7, 0, 7);
            _hPtypeKE->GetYaxis()->SetBinLabel(1, "mu");
            _hPtypeKE->GetYaxis()->SetBinLabel(2, "gamma");
            _hPtypeKE->GetYaxis()->SetBinLabel(3, "electron");
            _hPtypeKE->GetYaxis()->SetBinLabel(4, "neutron");
            _hPtypeKE->GetYaxis()->SetBinLabel(5, "proton");
            _hPtypeKE->GetYaxis()->SetBinLabel(6, "pion");
            _hPtypeKE->GetYaxis()->SetBinLabel(7, "kaon");

            _hSecondPtotVsPrimKE = tfdir.make<TH2D>("SecondPtotVsPrimKE",
                "Total momenta of secondaries vs Primary KE", 
                1000, 1E3, 1E8, 2000, 0., 2E6);
            _hShowerRadiusVsPrimKE = tfdir.make<TH2D>("primeKEvsShowerRadius", 
                "Primary KE vs shower radius", 1000, 1E3, 1E8, 300, 0., 300E3*1.42);
            _hNSecondariesVsPrimKE = tfdir.make<TH2D>("primeKEvsNSecondaries",
                "Primary KE vs number of secondaries", 1000, 1E3, 1E8, 200, 0., 200);
          }

          // And tree
          if (_saveTree) {
            makeTrees();
          }
        }


  CosmicCRY::~CosmicCRY() { }

  void CosmicCRY::generate( GenParticleCollection& genParts )
  {
    // Ref point, need to be here so that geometry info can be obtained
    // , but let's do this only once
    if (!_geomInfoObtained) {
      GeomHandle<Mu2eEnvelope> env;
      GeomHandle<WorldG4>  worldGeom;
      GeomHandle<ExtMonFNAL::ExtMon> extMonFNAL;
      GeomHandle<DetectorSystem> detsys;

      _envXmin = env->xmin();
      _envXmax = env->xmax();
      _envYmin = env->ymin();
      _envYmax = env->ymax();
      _envZmin = env->zmin();
      _envZmax = env->zmax();

      if (_refPointChoice == "TRACKER") 
        _cosmicReferencePointInMu2e = Hep3Vector(detsys->getOrigin().x(),
            _refY0, detsys->getOrigin().z());
      else if (_refPointChoice == "EXTMONFNAL") 
        _cosmicReferencePointInMu2e = 
          Hep3Vector(extMonFNAL->detectorCenterInMu2e().x(),
              _refY0, extMonFNAL->detectorCenterInMu2e().z());
      else if (_refPointChoice == "CALO") 
      {
        GeomHandle<Calorimeter> calorimeter;
        _cosmicReferencePointInMu2e = Hep3Vector(detsys->getOrigin().x(),
            _refY0, calorimeter->disk(0).geomInfo().origin().z());
      }
      else if (_refPointChoice == "UNDEFINED") 
        _cosmicReferencePointInMu2e = Hep3Vector(0., _refY0, 0.);
      
      _geomInfoObtained = true;
    }

    // std::cout << _cosmicReferencePointInMu2e << std::endl;
    // Getting CRY particles
    std::vector<CRYParticle*> *secondaries = new std::vector<CRYParticle*>;
    _cryGen->genEvent(secondaries);

    _pdgId0 = _cryGen->primaryParticle()->PDGid();
    _ke0 = _cryGen->primaryParticle()->ke(); // MeV
    _nSecondaries = secondaries->size();
    _t0 = _cryGen->timeSimulated();

    if (_doHistograms) {
      _hNSecondaries->Fill(secondaries->size());
      _hNSecondariesVsPrimKE->Fill(_ke0, _nSecondaries);
    }
    if (_saveTree) {
      _tCryPrimary->Fill();
    }

    double secondPtot = 0.;
    std::vector<CLHEP::Hep2Vector> secondXZ;

    for (unsigned j=0; j<secondaries->size(); j++) {
      CRYParticle* secondary = (*secondaries)[j];

      GlobalConstantsHandle<ParticleDataTable> pdt;
      const HepPDT::ParticleData& p_data = pdt->particle(secondary->PDGid()).ref();
      double mass = p_data.mass().value(); // in MeV

      double ke = secondary->ke(); // MeV by default in CRY
      double totalE = ke + mass;
      double totalP = safeSqrt(totalE * totalE - mass * mass);

      secondPtot += totalP;

      // Moving the CRY particle around: find all intersections with Mu2e
      // envelope, then set the position at the *first* intersection

      secondXZ.push_back(
          CLHEP::Hep2Vector(secondary->x() * 1000, secondary->y() * 1000));

      // Change coordinate system since y points upward, z points along
      // the beam line; which make cry(xyz) -> mu2e(zxy), uvw -> mu2e(zxy)
      CLHEP::Hep3Vector position(
          secondary->y() * 1000 + _cosmicReferencePointInMu2e.x(),
          secondary->z() * 1000 + _cosmicReferencePointInMu2e.y(),
          secondary->x() * 1000 + _cosmicReferencePointInMu2e.z()); // to mm
      CLHEP::HepLorentzVector mom4(totalP*secondary->v(), totalP*secondary->w(),
          totalP*secondary->u(), totalE);

      if (_projectToEnvelope) {
        _envIntersections.clear();
        calIntersections(position, CLHEP::Hep3Vector(secondary->v(),
              secondary->w(), secondary->u()));
      }

      if (_envIntersections.size() > 0) {
        int idx = 0;
        double highestY = _envIntersections.at(idx).y();
        for (unsigned i = 0; i < _envIntersections.size(); ++i) {
          if (_envIntersections.at(i).y() > highestY) {
            idx = i;
            highestY = _envIntersections.at(idx).y();
          }
        }

        genParts.push_back(
            GenParticle(static_cast<PDGCode::type>(secondary->PDGid()),
              GenId::cosmicCRY,
              _envIntersections.at(idx),
              mom4, secondary->t() - _t0));
      }

      if (_doHistograms) {
        _hXZ->Fill(position.x(), position.z());
        _hY->Fill(position.y());
        _hE->Fill(secondary->ke());
        // _hPlane->Fill(position.y());
        // Keep the original uvw order to have correct theta and phi.
        CLHEP::Hep3Vector momDir(secondary->u(), secondary->v(), secondary->w());
        _hTheta->Fill(momDir.theta());
        _hPhi->Fill(momDir.phi());

        _hPtot->Fill(totalP);
        _hTime->Fill(secondary->t() - _t0);
        _hLiveTime->Fill(secondary->t());
        _hPyOverPtot->Fill(secondary->w());
        switch (secondary->PDGid()) {
          case 13: // mu-
            _hNegMuKE->Fill(secondary->ke());
            _hPtypeKE->Fill(secondary->ke(), 0);
            break;
          case -13: // mu+
            _hPosMuKE->Fill(secondary->ke());
            _hPtypeKE->Fill(secondary->ke(), 0);
            break;
          case 22: // photon
            _hPtypeKE->Fill(secondary->ke(), 1);
            break;
          case -11: // e-
            _hPtypeKE->Fill(secondary->ke(), 2);
            break;
          case 11: // e+
            _hPtypeKE->Fill(secondary->ke(), 2);
            break;
          case 2112: // neutron
            _hPtypeKE->Fill(secondary->ke(), 3);
            break;
          case -2112: // neutron
            _hPtypeKE->Fill(secondary->ke(), 3);
            break;
          case 2212: // proton
            _hPtypeKE->Fill(secondary->ke(), 4);
            break;
          case -2212: // proton
            _hPtypeKE->Fill(secondary->ke(), 4);
            break;
          case 111: // pi0
            _hPtypeKE->Fill(secondary->ke(), 5);
            break;
          case 211: // pi+
            _hPtypeKE->Fill(secondary->ke(), 5);
            break;
          case -211: // pi-
            _hPtypeKE->Fill(secondary->ke(), 5);
            break;
          case 130: // k0 L
            _hPtypeKE->Fill(secondary->ke(), 6);
            break;
          case 310: // k0 S
            _hPtypeKE->Fill(secondary->ke(), 6);
            break;
          case 311: // k0
            _hPtypeKE->Fill(secondary->ke(), 6);
            break;
          case 321: // k+
            _hPtypeKE->Fill(secondary->ke(), 6);
            break;
          case -321: // k-
            _hPtypeKE->Fill(secondary->ke(), 6);
            break;
          default:
            _hPtypeKE->Fill(secondary->ke(), -1);
            break;
        }
      }

      if (_saveTree && j < _maxNSecondaries) {
        _pdgId1[j] = secondary->PDGid();
        _x1[j] = position.x();
        _y1[j] = position.y();
        _z1[j] = position.z();
        _t1[j] = secondary->t();
        _ke1[j] = secondary->ke();
        _px1[j] = totalP * secondary->v();
        _py1[j] = totalP * secondary->w();
        _pz1[j] = totalP * secondary->u();
        _ptot1[j] = totalP;

        // Keep the original uvw order to have correct theta and phi.
        CLHEP::Hep3Vector momDir(secondary->u(), secondary->v(), secondary->w()); 
        _theta1[j] = momDir.theta();
        _phi1[j] = momDir.phi();
      }

      if (_verbose > 1) {
        std::cout << "Secondary " << j 
          << ": " << CRYUtils::partName(secondary->id()) 
          << " (pdgId " << secondary->PDGid() << ")"
          << ", kinetic energy " << secondary->ke()  << " MeV"
          << ", position " 
          << "(" << secondary->x()
          << ", " << secondary->y()
          << ", " << secondary->z()
          << ") m"
          << ", time " << secondary->t() << " sec"
          << ", mass: " << mass
          << ", mom: " << totalP
          << ", mom dir.: " << secondary->u() <<", " << secondary->v()
          << ", " << secondary->w()
          << std::endl;

        // std::cout <<  genParts.back() << std::endl;
      }

      delete secondary;
    }

    if (_doHistograms) {
      _hSecondPtotVsPrimKE->Fill(_ke0, secondPtot);
      _hShowerRadiusVsPrimKE->Fill(_ke0, calMEC(secondXZ));
    }
    if (_saveTree) {
      _tCrySecondaries->Fill();
    }

    delete secondaries;
  }

  void CosmicCRY::makeTrees()
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("CosmicCRY");

    _tCryPrimary = tfdir.make<TTree>("cryPrimaryTree", "cryPrimaryTree");
    _tCryPrimary->Branch("evtId", &_evtId0);
    _tCryPrimary->Branch("pdgId", &_pdgId0);
    _tCryPrimary->Branch("ke", &_ke0);
    _tCryPrimary->Branch("nSecondaries", &_nSecondaries);

    _tCrySecondaries = tfdir.make<TTree>("crySecondariesTree",
        "crySecondariesTree");
    _tCrySecondaries->Branch("evtId", &_evtId0);
    _tCrySecondaries->Branch("ke0", &_ke0);
    _tCrySecondaries->Branch("nSecondaries", &_nSecondaries);
    _tCrySecondaries->Branch("pdgId", _pdgId1, "pdgId[nSecondaries]/I");
    _tCrySecondaries->Branch("x", _x1, "x[nSecondaries]/D");
    _tCrySecondaries->Branch("y", _y1, "y[nSecondaries]/D");
    _tCrySecondaries->Branch("z", _z1, "z[nSecondaries]/D");
    _tCrySecondaries->Branch("t", _t1, "t[nSecondaries]/D");
    _tCrySecondaries->Branch("ke", _ke1, "ke[nSecondaries]/D");
    _tCrySecondaries->Branch("px", _px1, "px[nSecondaries]/D");
    _tCrySecondaries->Branch("py", _py1, "py[nSecondaries]/D");
    _tCrySecondaries->Branch("pz", _pz1, "pz[nSecondaries]/D");
    _tCrySecondaries->Branch("ptot", _ptot1, "ptot[nSecondaries]/D");
    _tCrySecondaries->Branch("theta", _theta1, "theta[nSecondaries]/D");
    _tCrySecondaries->Branch("phi", _phi1, "phi[nSecondaries]/D");
  }

  void CosmicCRY::createSetupString()
  {
    if (_returnMuons) 
      _setupString.append("returnMuons 1 "); // must have the trailing white space
    else
      _setupString.append("returnMuons 0 ");

    if (_returnNeutrons) 
      _setupString.append("returnNeutrons 1 ");
    else
      _setupString.append("returnNeutrons 0 ");

    if (_returnProtons) 
      _setupString.append("returnProtons 1 ");
    else
      _setupString.append("returnProtons 0 ");

    if (_returnGammas) 
      _setupString.append("returnGammas 1 ");
    else
      _setupString.append("returnGammas 0 ");

    if (_returnElectrons) 
      _setupString.append("returnElectrons 1 ");
    else
      _setupString.append("returnElectrons 0 ");

    if (_returnPions) 
      _setupString.append("returnPions 1 ");
    else
      _setupString.append("returnPions 0 ");

    if (_returnKaons) 
      _setupString.append("returnKaons 1 ");
    else
      _setupString.append("returnKaons 0 ");

    char tmpStr[256];
    sprintf(tmpStr, "date %d-%d-%d ", _month, _day, _year);
    _setupString.append(tmpStr);

    sprintf(tmpStr, "latitude %f ", _latitude);
    _setupString.append(tmpStr);

    sprintf(tmpStr, "altitude %d ", _altitude);
    _setupString.append(tmpStr);

    sprintf(tmpStr, "subboxLength %f ", _subboxLength);
    _setupString.append(tmpStr);
  }

  double CosmicCRY::calMEC(std::vector<CLHEP::Hep2Vector> points)
  {
    double radius = 0.;
    unsigned int nPoints = points.size();

    if (nPoints <= 2)  {
      radius = 0.;
    }
    else
    {
      CLHEP::Hep2Vector center(-150E3, -150E3);
      double maxDistance = 0.;

      for (unsigned int i = 0; i < nPoints; ++i) {
        CLHEP::Hep2Vector diff = points.at(i) - center;
        if (diff.mag() > maxDistance) {
          maxDistance = diff.mag();
        }
      }

      radius = maxDistance;
    }
    return radius;
  }

  void CosmicCRY::calIntersections(CLHEP::Hep3Vector orig, CLHEP::Hep3Vector dir)
  {
    // roof: _envYmax, _envXmin, _envXmax, _envZmin, _envZmax
    // skip projection if the particle goes parallely to the plane
    if (dir.y() != 0.) {
      double t = (_envYmax - orig.y()) / dir.y();
      double x1 = dir.x() * t + orig.x();
      double z1 = dir.z() * t + orig.z();
      if (pointInBox(x1, z1, _envXmin, _envXmax, _envZmin, _envZmax)) {
        _envIntersections.push_back(CLHEP::Hep3Vector(x1, _envYmax, z1));
      }
    }

    //east: _envZmin, _envXmin, _envXmax, _envZmin, _envZmax
    if (dir.z() != 0.) {
      double t = (_envZmin - orig.z()) / dir.z();
      double x1 = dir.x() * t + orig.x();
      double y1 = dir.y() * t + orig.y();
      if (pointInBox(x1, y1, _envXmin, _envYmin, _envXmax, _envYmax)) {
        _envIntersections.push_back(CLHEP::Hep3Vector(x1, y1, _envZmin));
      }
    }

    //west: _envZmax, _envXmin, _envXmax, _envZmin, _envZmax
    if (dir.z() != 0.) {
      double t = (_envZmax - orig.z()) / dir.z();
      double x1 = dir.x() * t + orig.x();
      double y1 = dir.y() * t + orig.y();
      if (pointInBox(x1, y1, _envXmin, _envYmin, _envXmax, _envYmax)) {
        _envIntersections.push_back(CLHEP::Hep3Vector(x1, y1, _envZmax));
      }
    }

    //south: _envXmin, _envYmin, _envYmax, _envZmin, _envZmax
    if (dir.x() != 0.) {
      double t = (_envXmin - orig.x()) / dir.x();
      double z1 = dir.z() * t + orig.z();
      double y1 = dir.y() * t + orig.y();
      if (pointInBox(z1, y1, _envZmin, _envYmin, _envZmax, _envYmax)) {
        _envIntersections.push_back(CLHEP::Hep3Vector(_envXmin, y1, z1));
      }
    }

    //north: _envXmax, _envYmin, _envYmax, _envZmin, _envZmax
    if (dir.x() != 0.) {
      double t = (_envXmax - orig.x()) / dir.x();
      double z1 = dir.z() * t + orig.z();
      double y1 = dir.y() * t + orig.y();
      if (pointInBox(z1, y1, _envZmin, _envYmin, _envZmax, _envYmax)) {
        _envIntersections.push_back(CLHEP::Hep3Vector(_envXmax, y1, z1));
      }
    }
  }

  bool CosmicCRY::pointInBox(double x, double y, double x0, double y0,
      double x1, double y1)
  {
    bool ret = false;
    if ((x >= x0) && (x <= x1) && (y >= y0) && (y <= y1)) {
      ret = true;
    }
    return ret;
  }
}
