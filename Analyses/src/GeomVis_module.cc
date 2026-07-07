//
// Module to display geometry and event track
//
// Author: Zhengyun You, 03/21/2014
//

#include <iostream>
#include <iomanip>
#include <vector>

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
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
#include "TMath.h"
#include "TVector3.h"
#include "TGWindow.h"
#include "TRootEmbeddedCanvas.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoArb8.h"
#include "TGeoTube.h"

using namespace std;

namespace mu2e
{
  class GeomVis : public art::EDAnalyzer
  {
    typedef std::vector<art::InputTag> InputTags;
    InputTags stepInputs_;

    public:
      explicit GeomVis(fhicl::ParameterSet const &pset);
      virtual ~GeomVis() { }
      virtual void beginJob();
      void endJob();
      void analyze(const art::Event& e);


    private:
      void initGeom();
      void setVisible(TGeoNode* node);
      void drawLine(TVector3 pos1, TVector3 pos2, int pdg, int evt);
      TGeoManager* _geom;
      TGeoMedium*  _trkMat;
      int _trk;
      int _seg;

      fhicl::ParameterSet _pset;
      bool _first;

      TGMainFrame* _frame;
      TGHorizontalFrame *_mainFrame;
      TRootEmbeddedCanvas *_mainCanvas;
      TPad* _mainPad;
  };

  GeomVis::GeomVis(fhicl::ParameterSet const &pset)
    : art::EDAnalyzer(pset)
    , _pset(pset)
  {
    typedef std::vector<std::string> VS;
    const VS stepStrings(pset.get<VS>("stepInputs"));
    for(const auto& i : stepStrings) {
      stepInputs_.emplace_back(i);
    }
    _first = false;
    _trk = 0;
    _seg = 0;
  }

  void GeomVis::setVisible(TGeoNode* node)
  {
    TString nodeName = node->GetName();
    TString nameBase = TString(nodeName(0,4));
    if ( nameBase.Contains("TS") || nameBase.Contains("Coll") || nameBase.Contains("Pbar") ) {
      node->SetVisibility(1);
    }
    else {
      node->SetVisibility(0);
    }
  }

  void GeomVis::initGeom()
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

    for (int i1 = 0; i1 < hallAir->GetNdaughters(); i1++) {
      TGeoNode* node_lv1 = hallAir->GetNode(i1);
      setVisible(node_lv1);

      for (int i2 = 0; i2 < node_lv1->GetNdaughters(); i2++) {
        TGeoNode* node_lv2 = node_lv1->GetVolume()->GetNode(i2);
        setVisible(node_lv2);

        for (int i3 = 0; i3 < node_lv2->GetNdaughters(); i3++) {
          TGeoNode* node_lv3 = node_lv2->GetVolume()->GetNode(i3);
          setVisible(node_lv3);
        }
      }
    }
  }

  void GeomVis::drawLine(TVector3 pos1, TVector3 pos2, int pdg, int evt)
  {
    // pos 1, parent, -z, pos 2, self, +z, exchange if z's seqence not correct
    if (pos1.z() > pos2.z())
    {
      TVector3 tmp = pos1;
      pos1 = pos2;
      pos2 = tmp;
    }

    TVector3 center = 0.5*(pos1 + pos2);
    if (0) { cout<< "center " << endl;  center.Print(); }

    const TString base_name("evt_");
    TString name = base_name + Long_t(evt);
    name += TString("_trk") + Long_t(_trk);
    name += TString("_seg") + Long_t(_seg);
    _seg++;
    float dr = 1;

    int color = 1;
    if (pdg==PDGCode::anti_proton || pdg==PDGCode::pi_plus) color = 2;
    else if (pdg==PDGCode::pi_minus || pdg==PDGCode::mu_minus) color = 4;
    else if (pdg==PDGCode::proton) color = 3;
    color = _trk%10;

    TGeoVolume* top = _geom->GetTopVolume();
    int kline_shape = 2;
    if (kline_shape == 1)
    {
      float dz = 0.5*fabs(pos1.z() - pos2.z());
      TGeoArb8 *arb = new TGeoArb8(dz);
      arb->SetVertex(0, pos1.x()-center.x()-dr, pos1.y()-center.y()-dr);
      arb->SetVertex(1, pos1.x()-center.x()-dr, pos1.y()-center.y()+dr);
      arb->SetVertex(2, pos1.x()-center.x()+dr, pos1.y()-center.y()+dr);
      arb->SetVertex(3, pos1.x()-center.x()+dr, pos1.y()-center.y()-dr);
      arb->SetVertex(4, pos2.x()-center.x()-dr, pos2.y()-center.y()-dr);
      arb->SetVertex(5, pos2.x()-center.x()-dr, pos2.y()-center.y()+dr);
      arb->SetVertex(6, pos2.x()-center.x()+dr, pos2.y()-center.y()+dr);
      arb->SetVertex(7, pos2.x()-center.x()+dr, pos2.y()-center.y()-dr);

      TGeoVolume *vol = new TGeoVolume(name, arb, _trkMat);
      vol->SetLineColor(color);
      if (0) cout << "x " << center.x() << " y " << center.y() << endl;
      top->AddNode(vol, 1, new TGeoTranslation( center.x(), center.y(), center.z() ) );
    }
    else if (kline_shape == 2)
    {
      float dz = 0.5*(pos2 - pos1).Mag();
      TVector3 z_dir = pos2 - pos1;                                // direction of new z axis
      TVector3 x_dir = z_dir.Cross( TVector3(1.0, 2.0, 3.0) );     // a randomly chosen vector to get a new x axis perp to new z
      TVector3 y_dir = z_dir.Cross( x_dir );                       // cross to get new y

      TGeoTube *tube = new TGeoTube(0.0, dr, dz);

      TGeoRotation *rot = new TGeoRotation();
      float kdr = TMath::RadToDeg();
      rot->SetAngles(x_dir.Theta()*kdr, x_dir.Phi()*kdr, y_dir.Theta()*kdr, y_dir.Phi()*kdr, z_dir.Theta()*kdr, z_dir.Phi()*kdr);

      TGeoVolume *vol = new TGeoVolume(name, tube, 0);
      vol->SetLineColor(color);
      top->AddNode(vol, 1, new TGeoCombiTrans( center.x(), center.y(), center.z(), rot ) );
    }
  }

  void GeomVis::beginJob()
  {
    cout << "GeomVis beginJob " << endl;
    if (!gApplication) gApplication = new TApplication("GeomVis",0,0);
    initGeom();
  }

  void GeomVis::analyze(const art::Event& event)
  {
    int ks = 12;
    cout << "GeomVis analyze " << endl;

    if (_first) {
/*      int width = 1300, height = 700;
      _frame = new TGMainFrame(gClient->GetRoot(), width, height);
      _mainFrame = new TGHorizontalFrame(_frame, _frame->GetWidth(), _frame->GetHeight());
      _mainCanvas = new TRootEmbeddedCanvas("GeomVisCanvas",_mainFrame,width-200, height-100);
      _mainFrame->AddFrame(_mainCanvas, new TGLayoutHints(kLHintsTop));
      _mainFrame->MapWindow();

      _mainPad = new TPad("mainPad","Detector", 0, 0, 1, 1, 5,1,1);
      _mainPad->SetFillColor(1);
      _mainPad->Draw();
*/
      _first = false;
    }

    //TBrowser *_browser = new TBrowser();
    //gApplication->Run();

    _trk++;
    _seg = 0;
    for(const auto& tag : stepInputs_) {
      std::cout << tag << std::endl;
      auto ih = event.getValidHandle<StepPointMCCollection>(tag);
      for(const auto& i : *ih) {
        std::cout << "step (" << setw(ks) << i.position().x()
          << ", "  << setw(ks) << i.position().y()
          << ", "  << setw(ks) << i.position().z() << ")"
          << std::endl;

        double stepLength = i.stepLength();
        TVector3 posBegin(i.position().x(), i.position().y(), i.position().z());
        TVector3 momentum(i.momentum().x(), i.momentum().y(), i.momentum().z());
        momentum.SetMag(1.0);
        TVector3 posEnd = posBegin + momentum*stepLength;

        int pdg = PDGCode::anti_proton;
        int evt = event.id().event();
        drawLine(posBegin*0.1, posEnd*0.1, pdg, evt);
      }
    }
  }

  void GeomVis::endJob()
  {
    cout << "GeomVis endJob " << endl;

    //TBrowser *_browser = new TBrowser();

    //TH1D* h1 = new TH1D("h1", "h1", 100, 0, 100);
    //for (int i = 0; i < 100; i++) {
    //  h1->Fill(i, i);
    //}
    //h1->Draw();

    if (!_geom) {
      cout << "no geom found" << endl;
    }
    else {
      _geom->Export("mu2e_geom.root");
    }

    TGeoVolume* top = _geom->GetTopVolume();
    if (!top) cout << "top volume not found" << endl;
//    top->Draw("ogl");

    cout << "use button \"File->Quit ROOT\" to exit root" << endl;
//    gApplication->Run();
  }
}

using mu2e::GeomVis;
DEFINE_ART_MODULE(GeomVis)

