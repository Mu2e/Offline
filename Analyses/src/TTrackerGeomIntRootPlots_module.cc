//
// A plugin to do geometry plots using interactive root within the framework.
//
// $Id: TTrackerGeomIntRootPlots_module.cc,v 1.17 2013/10/21 20:44:04 genser Exp $
// $Author: genser $
// $Date: 2013/10/21 20:44:04 $
//
// Original author KLG based on Rob Kutschke's InteractiveRoot_plugin
//

#include "GeometryService/inc/GeomHandle.hh"
#include "TApplication.h"
#include "TArc.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyMarker.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"
#include "TTrackerGeom/inc/Support.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Plane.hh"
#include "TrackerGeom/inc/Panel.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib/pow.h"
#include "fhiclcpp/ParameterSet.h"
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

using cet::diff_of_squares;
using cet::sum_of_squares;

namespace mu2e {

  class Straw;

  class TTrackerGeomIntRootPlots : public art::EDAnalyzer {
  public:

    explicit TTrackerGeomIntRootPlots(fhicl::ParameterSet const& pset);
    virtual ~TTrackerGeomIntRootPlots() { }

//     virtual void analyze(Event const& ) = 0;
//     virtual void beginJob(){}
//     virtual void endJob(){}
//     virtual void endRun(Run const& ){}
//     virtual void beginSubRun(SubRun const& ){}
//     virtual void endSubRun(SubRun const& ){}
//     virtual void respondToOpenInputFile(FileBlock const& fb) {}
//     virtual void respondToCloseInputFile(FileBlock const& fb) {}
//     virtual void respondToOpenOutputFiles(FileBlock const& fb) {}
//     virtual void respondToCloseOutputFiles(FileBlock const& fb) {}


    void beginJob();
    void beginRun(art::Run const& );
    void endJob();

    // This is called for each event and is pure virtual
    void analyze(const art::Event& e);


    // other aux functions

  private:

    // Start: run time parameters

    // Pointers to histograms, ntuples, TGraphs.

    TH1F*         _hdummy;
    TCanvas*      _canvas;
    TApplication* _application;

    // aux funtions/data

    double _span;

    CLHEP::Hep3Vector _drawingOrigin;
    SimpleConfig const * _config;

    TTracker const *_ttracker;
    TDirectory const * _directory;
    char const * _dirname;
    TFile * _file;

    void drawEnvelopesSupport(bool dolabels);
    void drawPanel(bool dolabels);
    void drawPanelXZdetail(size_t dolabels,
                            double lx=-std::numeric_limits<double>::max(),
                            double ux= std::numeric_limits<double>::max(),
                            double lz=-std::numeric_limits<double>::max(),
                            double uz= std::numeric_limits<double>::max());
    void drawStraws(bool dolabels);
    void drawStrawsXZdetail(bool dolabels,
                            double lx=-std::numeric_limits<double>::max(),
                            double ux= std::numeric_limits<double>::max(),
                            double lz=-std::numeric_limits<double>::max(),
                            double uz= std::numeric_limits<double>::max());

    void drawManifolds(bool dolabels);

    double calculateDrawingAngle(double span, double radius);
    void printAndDraw(TLine* line, double x1, double y1, double x2, double y2);
    void drawArrowFromOrigin(double xt, double radius,
                             std::string const label,
                             short markerStyle=1, short markerColor=1, double xshift=5, double yshift=15);

    void labelPoint(double xt, double yt, std::string const label, double xshift=0, double yshift=0);

    void boundit(double& x,double lx, double ux);


  };

  TTrackerGeomIntRootPlots::TTrackerGeomIntRootPlots(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    // ROOT objects

    _hdummy(0),
    _canvas(0),
    _application(nullptr),
    _span(0),
    _drawingOrigin(CLHEP::Hep3Vector(0.,0.,0.)),
    _config(0),
    _ttracker(0),
    _directory(0),
    _dirname(""),
    _file(0){

  }

  void TTrackerGeomIntRootPlots::beginJob( ){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    gROOT->SetStyle("Plain");

    // Make a copy of the current style.
    TStyle* MyStyle = new TStyle(*gStyle);

    // Modify it.
    MyStyle->SetName("MyStyle");
    MyStyle->SetPadRightMargin(0.1);
    MyStyle->SetPadLeftMargin(0.1);
    MyStyle->SetPadBottomMargin(0.1);

    // This also resets the offsets.
    MyStyle->SetTitleXSize(0.025);
    MyStyle->SetTitleYSize(0.025);

    // These are scale factors that scale the default.
    MyStyle->SetTitleYOffset(1.9);
    MyStyle->SetTitleXOffset(1.15);

    //  MyStyle->SetTitleAlign(13);
    MyStyle->SetTitleX(0.15);

    MyStyle->SetTitleTextColor(kBlack);
    MyStyle->SetTitleFontSize(0.075);
    MyStyle->SetNdivisions(505,"X");
    MyStyle->SetNdivisions(505,"Y");
    MyStyle->SetLabelSize(0.035,"X");
    MyStyle->SetLabelSize(0.035,"Y");
    MyStyle->SetLabelOffset(0.01,"X");
    MyStyle->SetLabelOffset(0.005,"Y");

    MyStyle->SetOptStat("e");
    //MyStyle->SetOptStat("neuo");
    //MyStyle->SetOptStat("");
    MyStyle->SetStatFontSize(0.06);

    MyStyle->SetHistLineColor(kBlue);
    MyStyle->SetMarkerColor(kBlue);

    // Consult $ROOT_DIR/include/Gtypes.h
    MyStyle->SetMarkerStyle(kFullCircle);
    MyStyle->SetMarkerSize(0.5);

    MyStyle->cd(); //change the current style

    // Tell root to use current style for all objects
    gROOT->ForceStyle();

    // Create a root interactive environment.
    // This may be done only once per job.  If multiple modules need to make TCanvases, then
    // we need to find a way to coordinate this step; move it to the TFileService?
    int    tmp_argc(0);
    char** tmp_argv(0);
    _application = new TApplication( "noapplication", &tmp_argc, tmp_argv );

    // code to deal with a framework feature
    // first disable interactive behavior
    gROOT->SetBatch(1);
    static int const cpsize = 10;
    _canvas = tfs->make<TCanvas>("c0", "dummy canvas to force a directory creation",
                                 0, 0, cpsize, cpsize );
    _canvas->Close();
    delete _canvas;
    // back to interactive mode
    gROOT->SetBatch(0);
    std::cout << "We are in the directory named: " << gDirectory->GetName() << std::endl;
    _directory = gDirectory;
    _dirname = gDirectory->GetName();
    _file = gDirectory->GetFile();


  }


  void TTrackerGeomIntRootPlots::beginRun(art::Run const& ){

    // Get access to the master geometry system and its run time config.
    art::ServiceHandle<GeometryService> geom;
    _config = &(geom->config());

    GeomHandle<TTracker> ttracker;
    _ttracker = &*ttracker;

  }

  void TTrackerGeomIntRootPlots::analyze(const art::Event& event ) {
    // we do not do much here since we plot the detector geometry which does not change between events,
    // but it is a pure virtual function ...
  } // end analyze

  void TTrackerGeomIntRootPlots::endJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    std::cout << "We are in the directory named: " << gDirectory->GetName() << std::endl;
    // We need to make sure we are in the correct directory esp. before we write out the canvases
    _file->cd(_dirname);
    std::cout << "We are in the directory named: " << gDirectory->GetName() << std::endl;

    // Create a canvas for the xy view

    int cpsize = 1100;
    TString canvasName = "c1xy";
    _canvas = tfs->make<TCanvas>("c1", "xy", 0, 0, cpsize, cpsize );

    gDirectory->Append(_canvas);

    _span = 850.0; // mm

    TH1F* frame = _canvas->DrawFrame(
                                     _drawingOrigin.x(),
                                     _drawingOrigin.y(),
                                     _drawingOrigin.x()+_span,
                                     _drawingOrigin.y()+_span
                                     );


    // this is taken care by the style above
    //     frame->SetTitleSize(0.025,"X");
    //     frame->SetTitleSize(0.025,"Y");
    frame->SetTitle(";X(mm);Y(mm)");

    //if (false) {

    // here we do the geometry drawings

    drawEnvelopesSupport(true);
    drawPanel(true);
    drawStraws(true);
    drawManifolds(true);

    //}

    _canvas->Modified();
    _canvas->Update();

    std::cerr << "Double click on the Canvas to go to the next one" ;
    _canvas->WaitPrimitive();
    std::cerr << std::endl;

    _canvas->Print(canvasName+".png","png");


    // xz

    canvasName = "c2xz";
    _canvas = tfs->make<TCanvas>("c2", "xz", 0, 0, cpsize, cpsize );

    gDirectory->Append(_canvas);

    // xspan should be 161 for MECO 152 for 2 layer TTracker for now
    static double const xspan = 152; // mm
    static double const zspan = 10; // mm

    double _spanx = xspan;
    double _spanz = xspan;

    frame = _canvas->DrawFrame(
                               _drawingOrigin.x()-_spanx,
                               _drawingOrigin.z()-_spanz,
                               _drawingOrigin.x()+_spanx,
                               _drawingOrigin.z()+_spanz
                               );

    gPad->SetGridx(kTRUE);
    gPad->SetGridy(kTRUE);

    frame->SetTitle(";X(mm);Z(mm)");

    // here we do the geometry drawings

    drawPanelXZdetail(2);
    drawStrawsXZdetail(false);

    _canvas->Modified();
    _canvas->Update();

    std::cerr << "Double click on the Canvas to go to the next one" ;
    _canvas->WaitPrimitive();
    std::cerr << std::endl;
    _canvas->Print(canvasName+".png","png");


    // xz detail

    canvasName = "c3xzcenter";
    _canvas = tfs->make<TCanvas>("c3", "xz center detail", 0, 0, cpsize, cpsize );

    gDirectory->Append(_canvas);

    _spanx = zspan;
    _spanz = zspan;

    frame = _canvas->DrawFrame(
                               _drawingOrigin.x()-_spanx,
                               _drawingOrigin.z()-_spanz,
                               _drawingOrigin.x()+_spanx,
                               _drawingOrigin.z()+_spanz
                               );
    gPad->SetGridx(kTRUE);
    gPad->SetGridy(kTRUE);

    frame->SetTitle(";X(mm);Z(mm)");

    // here we do the geometry drawings

    drawPanelXZdetail(3,    -_spanx, _spanx, -_spanz, _spanz);
    drawStrawsXZdetail(true, -_spanx, _spanx, -_spanz, _spanz);

    std::cerr << "Double click on the Canvas to go to the next one" ;
    _canvas->WaitPrimitive();
    std::cerr << std::endl;
    _canvas->Print(canvasName+".png","png");

    // xz inner/left detail

    canvasName = "c4xzinner";
    _canvas = tfs->make<TCanvas>("c4", "xz inner detail", 0, 0, cpsize, cpsize );

    gDirectory->Append(_canvas);

    double _lx = -xspan;
    double _ux = -xspan + 2*zspan;
    double _lz = -zspan;
    double _uz =  zspan;

    frame = _canvas->DrawFrame(
                               _drawingOrigin.x()+_lx,
                               _drawingOrigin.z()+_lz,
                               _drawingOrigin.x()+_ux,
                               _drawingOrigin.z()+_uz
                               );
    gPad->SetGridx(kTRUE);
    gPad->SetGridy(kTRUE);

    frame->SetTitle(";X(mm);Z(mm)");

    // here we do the geometry drawings

    drawPanelXZdetail(4,   _lx, _ux, _lz, _uz);
    drawStrawsXZdetail(true,_lx, _ux, _lz, _uz);

    std::cerr << "Double click on the Canvas to close it" ;
    _canvas->WaitPrimitive();
    std::cerr << std::endl;
    _canvas->Print(canvasName+".png","png");

    // xz outer/right detail

    canvasName = "c5xzouter";
    _canvas = tfs->make<TCanvas>("c5", "xz outer detail", 0, 0, cpsize, cpsize );

    gDirectory->Append(_canvas);

    _lx =  xspan - 2*zspan;
    _ux =  xspan;
    _lz = -zspan;
    _uz =  zspan;

    frame = _canvas->DrawFrame(
                               _drawingOrigin.x()+_lx,
                               _drawingOrigin.z()+_lz,
                               _drawingOrigin.x()+_ux,
                               _drawingOrigin.z()+_uz
                               );
    gPad->SetGridx(kTRUE);
    gPad->SetGridy(kTRUE);

    frame->SetTitle(";X(mm);Z(mm)");

    // here we do the geometry drawings

    drawPanelXZdetail(5,   _lx, _ux, _lz, _uz);
    drawStrawsXZdetail(true,_lx, _ux, _lz, _uz);

    std::cerr << "Double click on the Canvas to close it" ;
    _canvas->WaitPrimitive();
    std::cerr << std::endl;
    _canvas->Print(canvasName+".png","png");

    gROOT->GetListOfCanvases()->Write();

    gROOT->GetListOfCanvases()->Delete();

  }

  void TTrackerGeomIntRootPlots::drawPanelXZdetail(size_t dolabels,
                                                    double lx, double ux,
                                                    double lz, double uz
                                                    ){

    TLine* line   = new TLine();
    line->SetLineColor(kRed);
    line->SetLineStyle(kDotted);

    // Draw the panel box in x,y

    const size_t iplane = 0;
    const Plane& plane = _ttracker->getPlane(iplane);

    const size_t isec = 0;
    const Panel& panel = plane.getPanel(isec);

    std::cout << "panel.boxHalfLengths : " <<
      panel.boxHalfLengths()[1] << " " <<
      panel.boxHalfLengths()[2] << " " <<
      panel.boxHalfLengths()[3] << " " <<
      panel.boxHalfLengths()[4] << " " << std::endl;

    std::cout << " panel.boxOffset : " <<
      panel.boxOffset().x() << " " <<
      panel.boxOffset().y() << " " <<
      panel.boxOffset().z() << " " << std::endl;

    std::cout << " _drawingOrigin : " <<
      _drawingOrigin.x() << " " <<
      _drawingOrigin.y() << " " <<
      _drawingOrigin.z() << " " << std::endl;

    // we do it arround "0" e.g. panel.boxOffset() as we have done the straws

    double x1 = _drawingOrigin.x()-panel.boxHalfLengths()[1];
    double x2 = x1;
    double z1 = _drawingOrigin.z()-panel.boxHalfLengths()[2];
    double z2 = _drawingOrigin.z()+panel.boxHalfLengths()[2];

    std::cout << "x1, z1 :" <<
      x1 << " " <<
      z1 << " " << std::endl;

    boundit (x1,lx,ux);
    boundit (z1,lz,uz);

    boundit (x2,lx,ux);
    boundit (z2,lz,uz);

    printAndDraw(line,x1,z1,x2,z2);

    double x3 = _drawingOrigin.x()+panel.boxHalfLengths()[1];
    double z3 = _drawingOrigin.z()+panel.boxHalfLengths()[2];

    boundit (x3,lx,ux);
    boundit (z3,lz,uz);

    printAndDraw(line,x2,z2,x3,z3);

    double x4 = _drawingOrigin.x()+panel.boxHalfLengths()[1];
    double z4 = _drawingOrigin.z()-panel.boxHalfLengths()[2];

    boundit (x4,lx,ux);
    boundit (z4,lz,uz);

    printAndDraw(line,x3,z3,x4,z4);

    printAndDraw(line,x1,z1,x4,z4);


    if (dolabels==2) {

      labelPoint(x1,z1,"PilE",0.,-7.);
      labelPoint(x2,z2,"PiuE",0., 3.);

      labelPoint(_drawingOrigin.x(),z1,"PilEm",-5.,-7.);
      labelPoint(_drawingOrigin.x(),z2,"PiuEm",-5., 3.);

      labelPoint(x3,z3,"PiuE",-10., 3.);
      labelPoint(x4,z4,"PilE",-10.,-7.);

    } else if (dolabels==3) {

      labelPoint(_drawingOrigin.x(),z1,"PilEm", -0.3,-0.5);
      labelPoint(_drawingOrigin.x(),z2,"PiuEm", -0.3, 0.25);

    } else if (dolabels==4) {

      labelPoint(x1,z1,"PilE",0.,0.);
      labelPoint(x2,z2,"PiuE",0.,0.);

    } else if (dolabels==5) {

      labelPoint(x1,z1,"PolE",0.,0.);
      labelPoint(x2,z2,"PouE",0.,0.);

    }

    delete line;

    return;

  }

  void TTrackerGeomIntRootPlots::drawStrawsXZdetail(bool dolabels,
                                                    double lx, double ux,
                                                    double lz, double uz
                                                    ){
    //    int _strawsPerManifold  = _config->getInt("ttracker.strawsPerManifold");

    TLine* line   = new TLine();
    line->SetLineStyle(kSolid);
    TArc* arc = new TArc();
    arc->SetFillStyle(0);
    arc->SetLineColor(kBlue);
    TPolyMarker* poly = new TPolyMarker();

    poly->SetMarkerSize(0.5);
    // poly->SetMarkerStyle(kFull);
    poly->SetMarkerColor(kOrange);

    const size_t iplane = 0;
    const Plane& plane = _ttracker->getPlane(iplane);
    const size_t isec = 0;
    const Panel& panel = plane.getPanel(isec);

    std::cout << "panel.boxOffset()  " <<
      panel.boxOffset() << std::endl;

    std::cout << "plane.origin()  " <<
      plane.origin() << std::endl;

    for ( int ilay=0; ilay<panel.nLayers(); ++ilay ) {

      const Layer& layer = panel.getLayer(ilay);

      if(ilay%2==0) {
        line->SetLineColor(kBlue);
      } else {
        line->SetLineColor(kMagenta);
      }

      // draw outer straw edges
      for ( int istr=0; istr<layer.nStraws(); ++istr ){

        const Straw& straw = layer.getStraw(istr);

        StrawDetail const& strawDetail = straw.getDetail();

        std::cout << "straw.getMidPoint() " << straw.getMidPoint() << std::endl;

        double sx = straw.getMidPoint().x() - panel.boxOffset().x();
        double sz = straw.getMidPoint().z() - panel.boxOffset().z();

        std::cout << "sx, sz :" <<
          sx << " " <<
          sz << " " << std::endl;

        if ( sx>ux || sz>uz || sx<lx || sz<lz ) continue;

        arc->DrawArc(sx,sz,strawDetail.outerRadius(),0.,360.,"only");
        poly->DrawPolyMarker( 1, &sx, &sz );

        if (dolabels) {

          std::ostringstream mlab("");

          mlab.width(6);
          mlab.fill(' ');
          mlab << "l" << ilay;
          mlab << "s" << istr;

          labelPoint(sx,sz,mlab.str(),-0.75,0.3);

        }

      }

    }

    delete line;
    delete arc;
    delete poly;

  }

  void TTrackerGeomIntRootPlots::drawStraws(bool dolabels) {

    int _strawsPerManifold = _config->getInt("ttracker.strawsPerManifold");

    TLine* line   = new TLine();
    line->SetLineStyle(kSolid);

    const size_t iplane = 0;
    const Plane& plane = _ttracker->getPlane(iplane);

    const size_t isec = 0;
    const Panel& panel = plane.getPanel(isec);

    for ( int ilay=0; ilay<panel.nLayers(); ++ilay ) {

      const Layer& layer = panel.getLayer(ilay);

      if(ilay%2==0) {
        line->SetLineColor(kBlue);
      } else {
        line->SetLineColor(kMagenta);
      }

      // draw outer straw edges
      for ( int istr=0; istr<layer.nStraws(); ++istr ){

        const Straw& straw = layer.getStraw(istr);

        StrawDetail const& strawDetail = straw.getDetail();

        double sx1 = _drawingOrigin.x() + straw.getMidPoint().x() + strawDetail.outerRadius();
        double sx2 = sx1;
        double sy1 = _drawingOrigin.y();
        double sy2 = _drawingOrigin.y() + strawDetail.halfLength();

        line->DrawLine( sx1, sy1, sx2, sy2 );

        // if (istr%_strawsPerManifold==0) {
        if (true) {
          // draw "lower" straw edge

          double sx1 = _drawingOrigin.x() + straw.getMidPoint().x() - strawDetail.outerRadius();
          double sx2 = sx1;
          double sy1 = _drawingOrigin.y();
          double sy2 = _drawingOrigin.y() + strawDetail.halfLength();

          line->SetLineColor(kGreen);
          line->DrawLine( sx1, sy1, sx2, sy2 );
          if(ilay%2==0) {
            line->SetLineColor(kBlue);
          } else {
            line->SetLineColor(kMagenta);
          }

          if (dolabels && istr%_strawsPerManifold==0 && istr!=0 && ilay==0) {

            // radius is sx2
            // drawArrowFromOrigin( sx2, sx2, "MT", kFullCircle, kOrange);

          }

        }

      }

    }
    delete line;
  }

  void TTrackerGeomIntRootPlots::drawManifolds(bool dolabels) {

    int _strawsPerManifold = _config->getInt("ttracker.strawsPerManifold");

    TLine* line   = new TLine();
    line->SetLineStyle(kSolid);

    const size_t iplane = 0;
    const Plane& plane = _ttracker->getPlane(iplane);

    const std::vector<double>& _manifoldHalfLengths = _ttracker->getManifoldHalfLengths();
    //const int _manifoldsPerEnd = _manifoldHalfLengths.size();

    const size_t isec = 0;
    const Panel& panel = plane.getPanel(isec);

    for ( int ilay=0; ilay<panel.nLayers(); ++ilay ) {

      const Layer& layer = panel.getLayer(ilay);

      if(ilay%2==0) {
        line->SetLineColor(kRed);
      } else {
        line->SetLineColor(kBlack);
      }

      // drawing manifolds (based on the first straw & manifold size)

      // deriving the strawGap

      static const int istr0=0;
      static const int istr1=1;
      Straw const& straw0 = layer.getStraw(istr0);
      Straw const& straw1 = layer.getStraw(istr1);
      StrawDetail const& strawDetail0 = straw0.getDetail();
      double _strawGap = straw1.getMidPoint().x() - straw0.getMidPoint().x() -
        strawDetail0.outerRadius()*2.0;

      double const dx = _manifoldHalfLengths.at(0) -
        strawDetail0.outerRadius()*_strawsPerManifold -
        _strawGap*(_strawsPerManifold-1)*0.5;

      std::cout << "_strawGap, dx : " <<
        _strawGap << ", " <<
        dx << std::endl;

      int iman = 0;
      for ( int istr=0; istr<layer.nStraws(); istr+=_strawsPerManifold ){

        ++iman;
        const Straw& straw = layer.getStraw(istr);

        StrawDetail const& strawDetail = straw.getDetail();

        double sx1 = _drawingOrigin.x() +
          straw.getMidPoint().x() - strawDetail.outerRadius() - dx;
        double sx2 = sx1;
        double sy1 = _drawingOrigin.y() + strawDetail.halfLength();
        double sy2 = sy1 + _manifoldHalfLengths.at(1)*2.0;

        std::cout << "sx1, sy1 :" <<
          sx1 << " " <<
          sy1 << " " << std::endl;

        printAndDraw(line, sx1, sy1, sx2, sy2 );

        double sx3 = sx1 + _manifoldHalfLengths.at(0)*2.0;

        printAndDraw(line,  sx1, sy1, sx3, sy1 );
        printAndDraw(line,  sx3, sy1, sx3, sy2 );
        printAndDraw(line,  sx1, sy2, sx3, sy2 );

        if (dolabels && ilay==0) {

          std::ostringstream mlab("");
          // mlab.str("");
          mlab.width(3);
          //mlab.fill('0');
          mlab << "m" << iman;

          labelPoint(sx1,sy1,mlab.str(),15.);

        }

      }

    }
    delete line;
  }


  void TTrackerGeomIntRootPlots::drawEnvelopesSupport(bool dolabels){

    TubsParams envelopeParams = _ttracker->getInnerTrackerEnvelopeParams();

    Support supportParams = _ttracker->getSupportParams();

    size_t const oldp = std::cout.precision();
    size_t const oldw = std::cout.width();
    std::cout << "Debugging tracker env envelopeParams ir,or,zhl,phi0,phimax:            " <<
      "   " <<
      std::fixed << std::setprecision(8) << std::setw(14) << envelopeParams.innerRadius() << ", " <<
      std::fixed << std::setprecision(8) << std::setw(14) << envelopeParams.outerRadius() << ", " <<
      std::fixed << std::setprecision(8) << std::setw(14) << envelopeParams.zHalfLength() << ", " <<
      std::fixed << std::setprecision(8) << std::setw(14) << envelopeParams.phi0()        << ", " <<
      std::fixed << std::setprecision(8) << std::setw(14) << envelopeParams.phiMax()
         << std::setprecision(oldp) << std::setw(oldw) << std::endl;

    TArc* arc = new TArc();
    arc->SetFillStyle(0);
    arc->SetLineColor(kCyan);


    double angle = calculateDrawingAngle(_span, envelopeParams.innerRadius());
    arc->DrawArc(_drawingOrigin.x(),_drawingOrigin.y(),
                 envelopeParams.innerRadius(),0.,angle,"only");

    arc->SetLineColor(kGreen);

    angle = calculateDrawingAngle(_span, supportParams.innerRadius());
    arc->DrawArc(_drawingOrigin.x(),_drawingOrigin.y(),
                 supportParams.innerRadius(),0.,angle,"only");

    angle = calculateDrawingAngle(_span, supportParams.outerRadius());
    arc->DrawArc(_drawingOrigin.x(),_drawingOrigin.y(),
                 supportParams.outerRadius(),0.,angle,"only");

//     double testCircleRadius = 720.0;
//     angle = calculateDrawingAngle(_span, testCircleRadius);
//     arc->DrawArc(_drawingOrigin.x(),_drawingOrigin.y(),
//                  testCircleRadius,0.,angle,"only");

    if (dolabels) {

      drawArrowFromOrigin( 75., envelopeParams.innerRadius(), "IER",kFullCircle,kRed);
      drawArrowFromOrigin(200., supportParams.innerRadius(),  "ISR",kFullCircle,kRed);
      drawArrowFromOrigin(300., supportParams.outerRadius(),  "OSR",kFullCircle,kRed);

    }

    return;

  }

  void TTrackerGeomIntRootPlots::drawPanel(bool dolabels){

    TLine* line   = new TLine();
    line->SetLineColor(kRed);
    line->SetLineStyle(kDotted);

    // Draw the panel box in x,y

    const size_t iplane = 0;
    const Plane& plane = _ttracker->getPlane(iplane);

    const size_t isec = 0;
    const Panel& panel = plane.getPanel(isec);

    std::cout << "panel.boxHalfLengths : " <<
      panel.boxHalfLengths()[1] << " " <<
      panel.boxHalfLengths()[2] << " " <<
      panel.boxHalfLengths()[3] << " " <<
      panel.boxHalfLengths()[4] << " " << std::endl;

    std::cout << " panel.boxOffset : " <<
      panel.boxOffset().x() << " " <<
      panel.boxOffset().y() << " " <<
      panel.boxOffset().z() << " " << std::endl;

    std::cout << " _drawingOrigin : " <<
      _drawingOrigin.x() << " " <<
      _drawingOrigin.y() << " " <<
      _drawingOrigin.z() << " " << std::endl;

    double x1 = -_drawingOrigin.x() + panel.boxOffset().x()-panel.boxHalfLengths()[1];
    double x2 = x1;
    double y1 = -_drawingOrigin.y();
    double y2 = -_drawingOrigin.y() + panel.boxHalfLengths()[4];

    std::cout << "x1, y1 :" <<
      x1 << " " <<
      y1 << " " << std::endl;

    printAndDraw(line,x1,y1,x2,y2);

    double x3 = -_drawingOrigin.x() + panel.boxOffset().x()+panel.boxHalfLengths()[1];
    double y3 = -_drawingOrigin.y() + panel.boxHalfLengths()[3];

    printAndDraw(line,x2,y2,x3,y3);

    double x4 = -_drawingOrigin.x() + panel.boxOffset().x()+panel.boxHalfLengths()[1];
    double y4 = -_drawingOrigin.y();

    printAndDraw(line,x3,y3,x4,y4);

    printAndDraw(line,x1,y1,x4,y4);


    if (dolabels) {

      double radius = sqrt(sum_of_squares(x2, y2));
      drawArrowFromOrigin( x2, radius, "PiE", kFullCircle, kRed);

      labelPoint( x2, y2, "#alpha",5.,-20.);

      drawArrowFromOrigin( x2, x2,     "PiM", kFullCircle, kRed,-50.);
      radius =  sqrt(sum_of_squares(x3, y3));
      drawArrowFromOrigin( x3, radius, "PoE", kFullCircle, kRed);
      drawArrowFromOrigin( x3, x3,     "PoM", kFullCircle, kRed);
    }

    delete line;

    return;

  }

  void TTrackerGeomIntRootPlots::labelPoint(double xt, double yt, std::string const label,
                                            double xshift, double yshift
                                            ){

    // Label the Point

    TLatex* text   = new TLatex();
    double ts = text->GetTextSize();
    //    std::cout << "Old Text Size : " << ts << std::endl;
    text->SetTextSize(ts*0.25);

    text->DrawLatex( xt+xshift, yt+yshift, label.c_str());

    delete text;

  }

  void TTrackerGeomIntRootPlots::printAndDraw(TLine* line, double x1, double y1, double x2, double y2) {

    std::cout << "x2, y2 :" <<
      x2 << " " <<
      y2 << " " << std::endl;

    line->DrawLine( x1, y1, x2, y2 );

  }

  void TTrackerGeomIntRootPlots::boundit(double& x,double lx, double ux) {

    if (x > ux) x = ux;
    if (x < lx) x = lx;

  }

  void TTrackerGeomIntRootPlots::drawArrowFromOrigin(double xt, double radius,
                                                     std::string const label,
                                                     short markerStyle, short markerColor,
                                                     double xshift, double yshift){

    // Label the Radii

    TText* text   = new TText();
    TPolyMarker* poly = new TPolyMarker();

    poly->SetMarkerStyle(markerStyle);
    poly->SetMarkerSize(0.75);
    poly->SetMarkerColor(markerColor);

    double ts = text->GetTextSize();
    //    std::cout << "Old Text Size : " << ts << std::endl;
    text->SetTextSize(ts*0.5);

    double yt = sqrt(diff_of_squares(radius, xt));

    text->DrawText( xt+xshift, yt+yshift, label.c_str());

    poly->DrawPolyMarker( 1, &xt, &yt   );

    TArrow* arrow = new TArrow();
    arrow->DrawArrow( _drawingOrigin.x(), _drawingOrigin.y(), xt, yt, 0.025, ">");
    //size & Option_t* option = "|>")

    delete poly;
    delete text;
    delete arrow;

  }

  double TTrackerGeomIntRootPlots::calculateDrawingAngle(double span, double radius) {

    if (span<0. || radius<0.) return 0.;
    double angle = (span<radius) ? 90.0 - acos(span/radius)*180./M_PI :
      90.0;

    return angle;

  }

}  // end namespace mu2e

using mu2e::TTrackerGeomIntRootPlots;
DEFINE_ART_MODULE(TTrackerGeomIntRootPlots);
