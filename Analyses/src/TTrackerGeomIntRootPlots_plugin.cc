//
// A plugin to do geometry plots using interactive root within the framework.
//
// $Id: TTrackerGeomIntRootPlots_plugin.cc,v 1.2 2010/10/01 15:44:35 genser Exp $
// $Author: genser $ 
// $Date: 2010/10/01 15:44:35 $
//
// Original author KLG based on Rob Kutschke's InteractiveRoot_plugin
//


// C++ includes.
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "GeneralUtilities/inc/pow.hh"

#include "GeometryService/inc/GeomHandle.hh"


#include "TrackerGeom/inc/Sector.hh"
#include "TrackerGeom/inc/Device.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "TTrackerGeom/inc/Support.hh"

// Mu2e includes.
//#include "ToyDP/inc/StepPointMCCollection.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1.h"

#include "TStyle.h"
#include "TFile.h"
#include "TROOT.h"

#include "TArc.h"
#include "TArrow.h"
#include "TText.h"
#include "TLine.h"

#include "TPolyMarker.h"

using namespace std;

namespace mu2e {

  class Straw;

  class TTrackerGeomIntRootPlots : public edm::EDAnalyzer {
  public:
    
    explicit TTrackerGeomIntRootPlots(edm::ParameterSet const& pset);
    virtual ~TTrackerGeomIntRootPlots() { }

//     virtual void analyze(Event const&, EventSetup const&) = 0;
//     virtual void beginJob(EventSetup const&){}
//     virtual void endJob(){}
//     virtual void endRun(Run const&, EventSetup const&){}
//     virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&){}
//     virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&){}
//     virtual void respondToOpenInputFile(FileBlock const& fb) {}
//     virtual void respondToCloseInputFile(FileBlock const& fb) {}
//     virtual void respondToOpenOutputFiles(FileBlock const& fb) {}
//     virtual void respondToCloseOutputFiles(FileBlock const& fb) {}


    void beginJob(edm::EventSetup const&);
    void beginRun(edm::Run const&, edm::EventSetup const&);
    void endJob();
 
    // This is called for each event and is pure virtual
    void analyze(const edm::Event& e, edm::EventSetup const&);


    // other aux functions

  private:

    // Start: run time parameters

    // Pointers to histograms, ntuples, TGraphs.

    TH1F*         _hdummy;
    TCanvas*      _canvas;
    TApplication* _application;
    
    // aux funtions/data

    double _span;

    CLHEP::Hep3Vector _mu2eDetectorOrigin;
    SimpleConfig const * _config;

    TTracker const *_ttracker;
    TDirectory const * _directory;
    char const * _dirname;
    TFile * _file;

    void setMu2eDetectorOrigin();

    void drawEnvelopesSupport(bool dolabels);
    void drawSector(bool dolabels);
    void drawStraws(bool dolabels);
    void drawStrawsXZ(bool dolabels);
    void drawManifolds(bool dolabels);

    double calculateDrawingAngle(double span, double radius);
    void printAndDraw(TLine* line, double x1, double y1, double x2, double y2);
    void drawArrowFromOrigin(double xt, double radius, string const label,
                   short markerStyle=1, short markerColor=1);

  };

  TTrackerGeomIntRootPlots::TTrackerGeomIntRootPlots(edm::ParameterSet const& pset) : 

    // ROOT objects
      
    _hdummy(0),
    _canvas(0),
    _application(0),
    _span(0),
    _mu2eDetectorOrigin(0),
    _config(0),
    _ttracker(0),
    _directory(0),
    _dirname(""),
    _file(0){

  }
  
  void TTrackerGeomIntRootPlots::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;

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
    cout << "We are in the directory named: " << gDirectory->GetName() << endl;
    _directory = gDirectory;
    _dirname = gDirectory->GetName();
    _file = gDirectory->GetFile();


  }


  void TTrackerGeomIntRootPlots::beginRun(edm::Run const&, edm::EventSetup const&){

    // Get access to the master geometry system and its run time config.
    edm::Service<GeometryService> geom;
    _config = &(geom->config());

    GeomHandle<TTracker> ttracker;
    _ttracker = &*ttracker;
    
  }

  void TTrackerGeomIntRootPlots::analyze(const edm::Event& event, edm::EventSetup const&) {
    // we do not do much here since we plot the detector geometry which does not change between events,
    // but it is a pure virtual function ...
  } // end analyze

  void TTrackerGeomIntRootPlots::endJob(){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;

    setMu2eDetectorOrigin();
    // we will overwrite it to be at "0"

    _mu2eDetectorOrigin = CLHEP::Hep3Vector(0.,0.,0.);

    std::cout << " _mu2eDetectorOrigin : " <<
      _mu2eDetectorOrigin.x() << " " <<
      _mu2eDetectorOrigin.y() << " " <<
      _mu2eDetectorOrigin.z() << " " << std::endl;


    cout << "We are in the directory named: " << gDirectory->GetName() << endl;
    // We need to make sure we are in the correct directory esp. before we write out the canvases
    _file->cd(_dirname);
    cout << "We are in the directory named: " << gDirectory->GetName() << endl;

    // Create a canvas for the xy view

    int cpsize = 1100;

    _canvas = tfs->make<TCanvas>("c1", "xy canvas", 0, 0, cpsize, cpsize );

    gDirectory->Append(_canvas);

    _span = 850.0; // mm

    TH1F* frame = _canvas->DrawFrame(
                                     _mu2eDetectorOrigin.x(),
                                     _mu2eDetectorOrigin.y(), 
                                     _mu2eDetectorOrigin.x()+_span, 
                                     _mu2eDetectorOrigin.y()+_span
                                     );

    // this is taken care by the style above
    //     frame->SetTitleSize(0.025,"X");
    //     frame->SetTitleSize(0.025,"Y");
    frame->SetTitle(";X(mm);Y(mm)");

    // here we do the geometry drawings

    drawEnvelopesSupport(true);    
    drawSector(true);
    drawStraws(true);
    drawManifolds(true);

    _canvas->Modified();
    _canvas->Update();

    cerr << "Double click on the Canvas to go to the next one" ;
    _canvas->WaitPrimitive();
    cerr << endl;

    // xz

    _canvas = tfs->make<TCanvas>("c2", "zy canvas, NOT FINISHED YET", 0, 0, cpsize, cpsize );

    gDirectory->Append(_canvas);

    double _spanx = 850.0; // mm
    double _spanz = 50.0;

    TubsParams envelopeParams = _ttracker->getTrackerEnvelopeParams();
    double xf = int(_mu2eDetectorOrigin.x()+envelopeParams.innerRadius/10.)*10.;

    std::cout << "xf :" << 
      xf << " " << std::endl;
    
    frame = _canvas->DrawFrame(
                               xf,
                               _mu2eDetectorOrigin.z(), 
                               _mu2eDetectorOrigin.x()+_spanx, 
                               _mu2eDetectorOrigin.z()+_spanz
                               );

    frame->SetTitle(";X(mm);Z(mm)");

    // here we do the geometry drawings

    drawStrawsXZ(true);

    cerr << "Double click on the Canvas to close it" ;
    _canvas->WaitPrimitive();
    cerr << endl;

    gROOT->GetListOfCanvases()->Write();
    gROOT->GetListOfCanvases()->Delete();

  }

  void TTrackerGeomIntRootPlots::drawEnvelopesSupport(bool dolabels){

    TubsParams envelopeParams = _ttracker->getTrackerEnvelopeParams();

    Support supportParams = _ttracker->getSupportParams();

    size_t const oldp = cout.precision();
    size_t const oldw = cout.width();
    cout << "Debugging tracker env envelopeParams ir,or,zhl,phi0,phimax:            " <<
      "   " << 
      fixed << setprecision(8) << setw(14) << envelopeParams.innerRadius << ", " <<
      fixed << setprecision(8) << setw(14) << envelopeParams.outerRadius << ", " <<
      fixed << setprecision(8) << setw(14) << envelopeParams.zHalfLength << ", " <<
      fixed << setprecision(8) << setw(14) << envelopeParams.phi0        << ", " <<
      fixed << setprecision(8) << setw(14) << envelopeParams.phiMax 
         << setprecision(oldp) << setw(oldw) << endl;

    TArc* arc = new TArc();
    arc->SetFillStyle(0);
    arc->SetLineColor(kCyan);


    double angle = calculateDrawingAngle(_span, envelopeParams.innerRadius);
    arc->DrawArc(_mu2eDetectorOrigin.x(),_mu2eDetectorOrigin.y(),
                 envelopeParams.innerRadius,0.,angle,"only");

    arc->SetLineColor(kGreen);

    angle = calculateDrawingAngle(_span, supportParams.innerRadius);
    arc->DrawArc(_mu2eDetectorOrigin.x(),_mu2eDetectorOrigin.y(),
                 supportParams.innerRadius,0.,angle,"only");

    angle = calculateDrawingAngle(_span, supportParams.outerRadius);
    arc->DrawArc(_mu2eDetectorOrigin.x(),_mu2eDetectorOrigin.y(),
                 supportParams.outerRadius,0.,angle,"only");

    if (dolabels) {

      drawArrowFromOrigin( 75., envelopeParams.innerRadius, "EiR",kFullCircle,kRed);
      drawArrowFromOrigin(200., supportParams.innerRadius,  "SiR",kFullCircle,kRed);
      drawArrowFromOrigin(300., supportParams.outerRadius,  "SoR",kFullCircle,kRed);
      
    }

    return;

  }

  void TTrackerGeomIntRootPlots::drawArrowFromOrigin(double xt, double radius, string const label, 
                                   short markerStyle, short markerColor){

    // Label the Radii

    TText* text   = new TText();
    TPolyMarker* poly = new TPolyMarker();

    poly->SetMarkerStyle(markerStyle);
    poly->SetMarkerSize(0.75);
    poly->SetMarkerColor(markerColor);
    
    double ts = text->GetTextSize();
    //    cout << "Old Text Size : " << ts << endl;
    text->SetTextSize(ts*0.5);

    double yt = sqrt(square(radius)-square(xt));

    text->DrawText( xt-5., yt+15., label.c_str());
    poly->DrawPolyMarker( 1, &xt, &yt   ); 

    TArrow* arrow = new TArrow();
    arrow->DrawArrow( _mu2eDetectorOrigin.x(), _mu2eDetectorOrigin.y(), xt, yt, 0.025, ">");
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


  void TTrackerGeomIntRootPlots::drawSector(bool dolabels){

    TLine* line   = new TLine();
    line->SetLineColor(kRed);
    line->SetLineStyle(kDotted);

    // Draw the sector box in x,y

    const size_t idev = 0;
    const Device& device = _ttracker->getDevice(idev);

    const size_t isec = 0;
    const Sector& sector = device.getSector(isec);

    std::cout << "sector.boxHalfLengths : " <<
      sector.boxHalfLengths()[1] << " " << 
      sector.boxHalfLengths()[2] << " " << 
      sector.boxHalfLengths()[3] << " " << 
      sector.boxHalfLengths()[4] << " " << std::endl;

    std::cout << " sector.boxOffset : " <<
      sector.boxOffset().x() << " " <<
      sector.boxOffset().y() << " " <<
      sector.boxOffset().z() << " " << std::endl;

    std::cout << " _mu2eDetectorOrigin : " <<
      _mu2eDetectorOrigin.x() << " " <<
      _mu2eDetectorOrigin.y() << " " <<
      _mu2eDetectorOrigin.z() << " " << std::endl;

    double x1 = -_mu2eDetectorOrigin.x() + sector.boxOffset().x()-sector.boxHalfLengths()[1];
    double x2 = x1;
    double y1 = -_mu2eDetectorOrigin.y();
    double y2 = -_mu2eDetectorOrigin.y() + sector.boxHalfLengths()[4];

    std::cout << "x1, y1 :" << 
      x1 << " " << 
      y1 << " " << std::endl;

    printAndDraw(line,x1,y1,x2,y2);

    double x3 = -_mu2eDetectorOrigin.x() + sector.boxOffset().x()+sector.boxHalfLengths()[1];
    double y3 = -_mu2eDetectorOrigin.y() + sector.boxHalfLengths()[3];

    printAndDraw(line,x2,y2,x3,y3);

    double x4 = -_mu2eDetectorOrigin.x() + sector.boxOffset().x()+sector.boxHalfLengths()[1];
    double y4 = -_mu2eDetectorOrigin.y();

    printAndDraw(line,x3,y3,x4,y4);

    printAndDraw(line,x1,y1,x4,y4);


    if (dolabels) {

      double radius = sqrt(square(x2)+square(y2));
      drawArrowFromOrigin( x2, radius, "PlE", kFullCircle, kRed);
      drawArrowFromOrigin( x2, x2,     "PlM", kFullCircle, kRed);
      radius =  sqrt(square(x3)+square(y3));
      drawArrowFromOrigin( x3, radius, "PsE", kFullCircle, kRed);
      drawArrowFromOrigin( x3, x3,     "PsM", kFullCircle, kRed);
    }

    delete line;

    return;

  }

  void TTrackerGeomIntRootPlots::drawStraws(bool dolabels) {

    int _strawsPerManifold  = _config->getInt("ttracker.strawsPerManifold");

    TLine* line   = new TLine();
    line->SetLineStyle(kSolid);

    const size_t idev = 0;
    const Device& device = _ttracker->getDevice(idev);

    const size_t isec = 0;
    const Sector& sector = device.getSector(isec);

    for ( int ilay=0; ilay<sector.nLayers(); ++ilay ) {

      const Layer& layer = sector.getLayer(ilay);
      
      if(ilay%2==0) {
        line->SetLineColor(kBlue);
      } else {
        line->SetLineColor(kMagenta);
      }
      
      // draw outer straw edges
      for ( int istr=0; istr<layer.nStraws(); ++istr ){

        const Straw& straw = layer.getStraw(istr);

        StrawDetail const& strawDetail = straw.getDetail();

        double sx1 = _mu2eDetectorOrigin.x() + straw.getMidPoint().x() + strawDetail.outerRadius();
        double sx2 = sx1;
        double sy1 = _mu2eDetectorOrigin.y();
        double sy2 = _mu2eDetectorOrigin.y() + strawDetail.halfLength();

        line->DrawLine( sx1, sy1, sx2, sy2 );

        // if (istr%_strawsPerManifold==0) {
        if (true) {
          // draw "lower" straw edge
          
          double sx1 = _mu2eDetectorOrigin.x() + straw.getMidPoint().x() - strawDetail.outerRadius();
          double sx2 = sx1;
          double sy1 = _mu2eDetectorOrigin.y();
          double sy2 = _mu2eDetectorOrigin.y() + strawDetail.halfLength();

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

  void TTrackerGeomIntRootPlots::drawStrawsXZ(bool dolabels) {

    int _strawsPerManifold  = _config->getInt("ttracker.strawsPerManifold");

    TLine* line   = new TLine();
    line->SetLineStyle(kSolid);

    TArc* arc = new TArc();
    arc->SetFillStyle(0);
    arc->SetLineColor(kBlue);

    const size_t idev = 0;
    const Device& device = _ttracker->getDevice(idev);

    const size_t isec = 0;
    const Sector& sector = device.getSector(isec);

    for ( int ilay=0; ilay<sector.nLayers(); ++ilay ) {

      const Layer& layer = sector.getLayer(ilay);
      
      if(ilay%2==0) {
        line->SetLineColor(kBlue);
      } else {
        line->SetLineColor(kMagenta);
      }
      
      // draw outer straw edges
      for ( int istr=0; istr<layer.nStraws(); ++istr ){

        const Straw& straw = layer.getStraw(istr);

        StrawDetail const& strawDetail = straw.getDetail();

        double sx = _mu2eDetectorOrigin.x() + straw.getMidPoint().x();
        double sz = _mu2eDetectorOrigin.z() + straw.getMidPoint().z();

        std::cout << "sx, sz :" << 
          sx << " " << 
          sz << " " << std::endl; 


        //xf :370.00000000
        //sx, sz :382.50000000 -1460.00000000


        arc->DrawArc(sx,sz,strawDetail.outerRadius(),0.,360.,"only");

        if (dolabels && istr%_strawsPerManifold==0 && istr!=0 && ilay==0) {
            

        }
 
      }

    }

    delete line;
  }

  void TTrackerGeomIntRootPlots::drawManifolds(bool dolabels) {

    int _strawsPerManifold = _config->getInt("ttracker.strawsPerManifold");

    TLine* line   = new TLine();
    line->SetLineStyle(kSolid);

    const size_t idev = 0;
    const Device& device = _ttracker->getDevice(idev);

    const std::vector<double>& _manifoldHalfLengths = _ttracker->getManifoldHalfLengths();
    //const int _manifoldsPerEnd = _manifoldHalfLengths.size();

    const size_t isec = 0;
    const Sector& sector = device.getSector(isec);

    for ( int ilay=0; ilay<sector.nLayers(); ++ilay ) {

      const Layer& layer = sector.getLayer(ilay);
      
      if(ilay%2==0) {
        line->SetLineColor(kRed);
      } else {
        line->SetLineColor(kBlack);
      }
      
      // draw manifolds (based on the first straw & manifold size)            ----------------

      // get strawGap

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

        double sx1 = _mu2eDetectorOrigin.x() + 
          straw.getMidPoint().x() - strawDetail.outerRadius() - dx;
        double sx2 = sx1;
        double sy1 = _mu2eDetectorOrigin.y() + strawDetail.halfLength();
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
        
          ostringstream mlab("");
          // mlab.str("");
          mlab.width(3);
          //mlab.fill('0');
          mlab << "m" << iman;
          
          double radius = sqrt(square(sx1)+square(sy1));
          drawArrowFromOrigin( sx1, radius, mlab.str(), kFullCircle, kOrange);
        }

      }

    }
    delete line;
  }


  void TTrackerGeomIntRootPlots::printAndDraw(TLine* line, double x1, double y1, double x2, double y2) {

    std::cout << "x2, y2 :" << 
      x2 << " " << 
      y2 << " " << std::endl;
    
    line->DrawLine( x1, y1, x2, y2 );
    
  }

  void TTrackerGeomIntRootPlots::setMu2eDetectorOrigin() {

    // Dimensions of the world.
    vector<double> worldHLen;
    _config->getVectorDouble("world.halfLengths", worldHLen, 3);

    // Floor thickness.
    double floorThick = _config->getDouble("hall.floorThick");

    // Top of the floor in G4 world coordinates.
    double yFloor = -worldHLen[1] + floorThick;

    // The height above the floor of the y origin of the Mu2e coordinate system.
    double yOriginHeight = _config->getDouble("world.mu2eOrigin.height" )*CLHEP::mm;

    // Position of the origin of the mu2e coordinate system
    CLHEP::Hep3Vector _mu2eOrigin = 
      CLHEP::Hep3Vector( 
                        _config->getDouble("world.mu2eOrigin.xoffset")*CLHEP::mm,
                        yFloor + yOriginHeight,
                        _config->getDouble("world.mu2eOrigin.zoffset")*CLHEP::mm
                        );
    
    // Origin used to construct the MECO detector.
    // Magic number to fix:

    _mu2eDetectorOrigin = _mu2eOrigin + CLHEP::Hep3Vector( -3904., 0., 12000.);
  
    std::cout << " _mu2eOrigin : " <<
      _mu2eOrigin.x() << " " <<
      _mu2eOrigin.y() << " " <<
      _mu2eOrigin.z() << " " << std::endl;

    std::cout << " _mu2eDetectorOrigin : " <<
      _mu2eDetectorOrigin.x() << " " <<
      _mu2eDetectorOrigin.y() << " " <<
      _mu2eDetectorOrigin.z() << " " << std::endl;

  }

}  // end namespace mu2e

using mu2e::TTrackerGeomIntRootPlots;
DEFINE_FWK_MODULE(TTrackerGeomIntRootPlots);
