//
// Module to analyze tracks in TS, for pbar and muon separation
//
// Author: Zhengyun You, 03/27/2014
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "TApplication.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGWindow.h"
#include "TRootEmbeddedCanvas.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoArb8.h"
#include "TGeoTube.h"
#include "TMath.h"

using namespace std;

namespace mu2e
{
  class TSTrackAna : public art::EDAnalyzer
  {
    typedef std::vector<art::InputTag> InputTags;
    InputTags stepInputs_;
    typedef SimParticleCollection::key_type key_type;

    public:
      explicit TSTrackAna(fhicl::ParameterSet const &pset);
      virtual ~TSTrackAna() { }
      virtual void beginJob();
      void endJob();
      void analyze(const art::Event& e);

    private:
      void initGeom();
      double getTheta(const TVector3 &pos);
      void addIntersection(int ts, const StepPointMC &step);
      int  geth2Id(int ts, int i);
      void fill2DHist(int ts);

      TGeoManager* _geom;
      TGeoMedium*  _trkMat;
      int _trk;
      int _seg;

      TVector3 _originTS2;
      double _TS2R;
      double _deltaMax;
      double _cutThetaBegin;
      double _cutThetaEnd;
      double _cutThetaStep;
      double _cutTheta;
      double _cutZBegin;
      double _cutZEnd;
      double _cutZStep;
      double _cutZ;
      double _cutXBegin;
      double _cutXEnd;
      double _cutXStep;
      double _cutX;
      int _verbosity;

      fhicl::ParameterSet _pset;
      bool _first;
      int _ks;

      std::map<int, TVector3> _mapTrkIntersection;
      std::map<int, double>   _mapTrkDelta;
      std::map<int, TH2F*> _maph2YvsR;
  };

  TSTrackAna::TSTrackAna(fhicl::ParameterSet const &pset)
    : art::EDAnalyzer(pset)
    , _pset(pset)
  {
    typedef std::vector<std::string> VS;
    const VS stepStrings(pset.get<VS>("stepInputs"));
    for(const auto& i : stepStrings) {
      stepInputs_.emplace_back(i);
    }
    _cutThetaBegin = pset.get<double>("cutThetaBegin", 0.0);
    _cutThetaEnd   = pset.get<double>("cutThetaEnd",   90.0);
    _cutThetaStep  = pset.get<double>("cutThetaStep",  1.0);
    _cutZBegin = pset.get<double>("cutZBegin", -4044.0);
    _cutZEnd   = pset.get<double>("cutZEnd",   -2929.0);
    _cutZStep  = pset.get<double>("cutZStep",  5.0);
    _cutXBegin = pset.get<double>("cutXBegin", 975.0);
    _cutXEnd   = pset.get<double>("cutXEnd",   0.0);
    _cutXStep  = pset.get<double>("cutXStep",  -5.0);
    _verbosity = pset.get<int>("verbosity", 0);

    _first = false;
    _trk = 0;
    _seg = 0;
    _ks = 12;
    _originTS2 = TVector3(975.0, 0.0, -2929.0);
    _TS2R = 2929.0;
    _deltaMax = 50.0; // mm, distance from step point to cut plane

    art::ServiceHandle<art::TFileService> tfs;
    int nbin = 100;
    double rVac = 250.0;
    TString h2YvsRBaseName("h2YvsR_");
    TString h2YvsRBaseTitle("Y vs R at cut plane ");

    int nZ = int((_cutZEnd-_cutZBegin)/_cutZStep);
    for (int iZ = 0; iZ <= nZ; iZ++) {
      double z = _cutZBegin + _cutZStep*iZ;
      TString h2YvsRName = h2YvsRBaseName;
      h2YvsRName += geth2Id(1, z);
      TString h2YvsRTitle = h2YvsRBaseTitle;
      h2YvsRTitle += " z=";
      h2YvsRTitle += int(z);
      TH2F* h2YvsR = tfs->make<TH2F>( h2YvsRName, h2YvsRTitle, nbin, -_originTS2.Z()-rVac, -_originTS2.Z()+rVac, nbin, -rVac, rVac);
      _maph2YvsR.insert(std::pair<int, TH2F*>(geth2Id(1, int(z)), h2YvsR));
    }

    int nTheta = int((_cutThetaEnd-_cutThetaBegin)/_cutThetaStep);
    for (int iTheta = 0; iTheta <= nTheta; iTheta++) {
      double theta = _cutThetaBegin + _cutThetaStep*iTheta;
      TString h2YvsRName = h2YvsRBaseName;
      h2YvsRName += geth2Id(2, int(theta));
      TString h2YvsRTitle = h2YvsRBaseTitle;
      h2YvsRTitle += " theta=";
      h2YvsRTitle += int(theta);
      TH2F* h2YvsR = tfs->make<TH2F>( h2YvsRName, h2YvsRTitle, nbin, -_originTS2.Z()-rVac, -_originTS2.Z()+rVac, nbin, -rVac, rVac);
      _maph2YvsR.insert(std::pair<int, TH2F*>(geth2Id(2, int(theta)), h2YvsR));
    }

    int nX = int((_cutXEnd-_cutXBegin)/_cutXStep);
    for (int iX = 0; iX <= nX; iX++) {
      double x = _cutXBegin + _cutXStep*iX;
      TString h2YvsRName = h2YvsRBaseName;
      h2YvsRName += geth2Id(3, x);
      TString h2YvsRTitle = h2YvsRBaseTitle;
      h2YvsRTitle += " x=";
      h2YvsRTitle += int(x);
      TH2F* h2YvsR = tfs->make<TH2F>( h2YvsRName, h2YvsRTitle, nbin, -_originTS2.Z()-rVac, -_originTS2.Z()+rVac, nbin, -rVac, rVac);
      _maph2YvsR.insert(std::pair<int, TH2F*>(geth2Id(3, int(x)), h2YvsR));
    }

  }

  void TSTrackAna::initGeom()
  {
    _geom = TGeoManager::Import("mu2e_common.gdml");
    if (!_geom) cout << "geometry gdml file does not exist" << endl;
    //_geom->Export("mu2e.C");

    TGeoMaterial* vacuum = new TGeoMaterial("vacuum",200,100,20);
    _trkMat = new TGeoMedium("TrkMat",1,vacuum);

    TGeoVolume* top = _geom->GetTopVolume();
    TGeoVolume *hallAir = 0;
    for (int i0 = 0; i0 < top->GetNdaughters(); i0++) {
      TGeoNode* node_lv0 = top->GetNode(i0);
      if ( TString(node_lv0->GetName()).Contains("HallAir")) {
        hallAir = node_lv0->GetVolume();
        cout << "hallAir name " << hallAir->GetName() << endl;
      }
      else {
        node_lv0->SetVisibility(0);
      }
    }
  }

  void TSTrackAna::beginJob()
  {
    cout << "TSTrackAna beginJob " << endl;
    //if (!gApplication) gApplication = new TApplication("TSTrackAna",0,0);
    //initGeom();
  }

  int TSTrackAna::geth2Id(int ts, int i)
  {
    return ts*10000+i;
  }

  double TSTrackAna::getTheta(const TVector3 &pos)
  {
    TVector3 localPos = pos-_originTS2;
    localPos.SetY(0.0);
    double theta = localPos.Theta();
    if (localPos.x() < 0) theta *= -1;
    double r = localPos.Mag();
    if (_verbosity >= 3) {
      std::cout << "(" << setw(_ks) << localPos.X() << ", " << setw(_ks) << localPos.Z() << ") "
        << " r=" << setw(_ks) << r << " theta=" << setw(_ks) << theta << std::endl;
    }
    return theta;
  }

  void TSTrackAna::addIntersection(int ts, const StepPointMC &step)
  {
    TVector3 pos(step.position().x(), step.position().y(), step.position().z());
    double theta = getTheta(pos);
    double z = pos.Z();
    double x = pos.X();

    double delta = 9e99;
    if (ts == 1)      delta = _cutZ - z;
    else if (ts == 2) delta = -1.0*(_cutTheta*TMath::DegToRad() - theta)*_TS2R;
    else if (ts == 3) delta = -1.0*(_cutX - x);

    if ((ts == 1 || ts ==3) && fabs(delta) > _deltaMax) return;
    if (ts == 2 && (theta > TMath::PiOver2()+0.01 || theta < -0.01 || fabs(delta) > _deltaMax)) return;

    TVector3 dir(step.momentum().x(), step.momentum().y(), step.momentum().z());
    dir.SetMag(1.0);
    double distance = 0;
    if (ts == 1)      distance = delta/fabs(dir.Z());
    else if (ts == 2) distance = delta/sqrt(dir.X()*dir.X()+dir.Z()*dir.Z());
    else if (ts == 3) distance = delta/fabs(dir.X());

    TVector3 cutPos = pos + distance*dir;
    if (_verbosity >= 1) {
      cout << "pos (" << setw(_ks) << pos.X() << ", "  << setw(_ks) << pos.Y() << ", "  << setw(_ks) << pos.Z() << ")" << std::endl;
      dir.Print();
      cutPos.Print();
    }

    key_type trackId = step.trackId();
    int id = trackId.asInt();

    if (_mapTrkDelta.find(id) == _mapTrkDelta.end()) {
      _mapTrkIntersection.insert(std::pair<int, TVector3>(id, cutPos));
      _mapTrkDelta.insert(std::pair<int, double>(id, fabs(delta)));
      if (_verbosity >= 1) cout << "add new for id " << id << "with delta " << delta << endl;
    }
    else {
      if (fabs(delta) < _mapTrkDelta[id]) {
        _mapTrkIntersection[id] = cutPos;
        _mapTrkDelta[id] = fabs(delta);
        if (_verbosity >= 1) cout << "replace id " << id << " with delta " << delta << endl;
      }
    }
  }

  void TSTrackAna::fill2DHist(int ts)
  {
    int h2Id = 0;
    if (ts == 1)      h2Id = geth2Id(ts, int(_cutZ));
    else if (ts == 2) h2Id = geth2Id(ts, int(_cutTheta));
    else if (ts == 3) h2Id = geth2Id(ts, int(_cutX));

    TH2F* h2YvsR = _maph2YvsR[h2Id];
    if (!h2YvsR) cout << "h2id " << h2Id << " does not exist in _maph2YvsR" << endl;

    for(const auto& i: _mapTrkIntersection) {
      TVector3 localPos = i.second-_originTS2;
      localPos.SetY(0);
      double r = _originTS2.Z();
      if (ts == 1)      r = localPos.X();
      else if (ts == 2) r = localPos.Mag();
      else if (ts == 3) r = localPos.Z();
      double y = i.second.Y();
      if (_verbosity >= 1) cout << "r=" << setw(_ks) << r << " y=" << setw(_ks) << y << endl;
      h2YvsR->Fill(r, y);
    }
    _mapTrkIntersection.clear();
    _mapTrkDelta.clear();
  }

  void TSTrackAna::analyze(const art::Event& event)
  {
    if (_first) {
      _first = false;
    }

    _trk++;
    _seg = 0;

    int nZ = int((_cutZEnd-_cutZBegin)/_cutZStep);
    for (int iZ = 0; iZ <= nZ; iZ++) {
      _cutZ = _cutZBegin + _cutZStep*iZ;

      for(const auto& tag : stepInputs_) {
        //std::cout << tag << std::endl;
        auto ih = event.getValidHandle<StepPointMCCollection>(tag);
        for(const auto& i : *ih) {
          if (_verbosity >= 2) {
            std::cout << "Z step (" << setw(_ks) << i.position().x()
              << ", "  << setw(_ks) << i.position().y()
              << ", "  << setw(_ks) << i.position().z() << ")"
              << std::endl;
          }
          addIntersection(1, i);
        }
      }
      fill2DHist(1);  //TS1
    }

    int nTheta = int((_cutThetaEnd-_cutThetaBegin)/_cutThetaStep);
    for (int iTheta = 0; iTheta <= nTheta; iTheta++) {
      _cutTheta = _cutThetaBegin + _cutThetaStep*iTheta;

      for(const auto& tag : stepInputs_) {
        //std::cout << tag << std::endl;
        auto ih = event.getValidHandle<StepPointMCCollection>(tag);
        for(const auto& i : *ih) {
          if (_verbosity >= 2) {
            std::cout << "Theta step (" << setw(_ks) << i.position().x()
              << ", "  << setw(_ks) << i.position().y()
              << ", "  << setw(_ks) << i.position().z() << ")"
              << std::endl;
          }
          addIntersection(2, i);
        }
      }
      fill2DHist(2);  //TS2
    }

    int nX = int((_cutXEnd-_cutXBegin)/_cutXStep);
    for (int iX = 0; iX <= nX; iX++) {
      _cutX = _cutXBegin + _cutXStep*iX;

      for(const auto& tag : stepInputs_) {
        //std::cout << tag << std::endl;
        auto ih = event.getValidHandle<StepPointMCCollection>(tag);
        for(const auto& i : *ih) {
          if (_verbosity >= 2) {
            std::cout << "X step (" << setw(_ks) << i.position().x()
              << ", "  << setw(_ks) << i.position().y()
              << ", "  << setw(_ks) << i.position().z() << ")"
              << std::endl;
          }
          addIntersection(3, i);
        }
      }
      fill2DHist(3);  //TS3
    }

  }

  void TSTrackAna::endJob()
  {
    cout << "TSTrackAna endJob " << endl;

    if (!_geom) {
      cout << "no geom found" << endl;
    }
    else {
      _geom->Export("mu2e_geom.root");

      TGeoVolume* top = _geom->GetTopVolume();
      if (!top) cout << "top volume not found" << endl;
//    top->Draw("ogl");
    }
  }
}

using mu2e::TSTrackAna;
DEFINE_ART_MODULE(TSTrackAna)

