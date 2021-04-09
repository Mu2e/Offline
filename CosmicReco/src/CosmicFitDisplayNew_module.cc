//TGeo:
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoMatrix.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

//Mu2e Data Prods:
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"

#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/ComboHit.hh"

#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
// Mu2e Utilities
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
#include "Mu2eUtilities/inc/ParametricFit.hh"

//Mu2e Tracker Geom:
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"

// Mu2e diagnostics
#include "TrkDiag/inc/ComboHitInfo.hh"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"

//Cosmics:
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// ROOT incldues
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TH2D.h"
#include "TF1.h"
#include "TH3D.h"
#include "Rtypes.h"
#include "TApplication.h"
#include "TArc.h"
#include "TTUBE.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TText.h"
#include "TRotMatrix.h"


using namespace std; 

namespace mu2e 
{
  class CosmicFitDisplayNew : public art::EDAnalyzer {
    public:
	    struct Config{
	          using Name=fhicl::Name;
	          using Comment=fhicl::Comment;
	          fhicl::Atom<bool> mcdiag{Name("mcdiag"), Comment("set on for MC info"),true};
	          fhicl::Atom<art::InputTag> chtag{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
	          fhicl::Atom<art::InputTag> tctag{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
	          fhicl::Atom<art::InputTag> sttag{Name("CosmicTrackSeedCollection"),Comment("tag for cosmci track seed collection")};
	          fhicl::Atom<bool> doDisplay{Name("doDisplay"),Comment("use display"), false};
	          fhicl::Atom<bool> clickToAdvance{Name("clickToAdvance"),Comment("next event"), false};
	     };
        typedef art::EDAnalyzer::Table<Config> Parameters;
        explicit CosmicFitDisplayNew(const Parameters& conf);
        virtual ~CosmicFitDisplayNew();
        virtual void beginJob() override;
        virtual void analyze(const art::Event& e) override;
    private: 
      Config _conf;
      bool _mcdiag;
      Int_t _evt; 

      // The module label of this instance of this module.
      std::string moduleLabel_;
      const ComboHitCollection* _chcol;
      const CosmicTrackSeedCollection* _coscol;
      //For Event Displays:
      TApplication* application_;
      TDirectory*   directory_ = nullptr;
      TCanvas*      canvas_ = nullptr;

      art::InputTag   _chtag;//combo
      art::InputTag   _tctag;//timeclusters
      art::InputTag   _sttag;//Straight tracks
      bool doDisplay_;
      bool clickToAdvance_;

      void plot2d(const art::Event& evt);
      void plot3dXYZ(const art::Event& evt);
      std::vector<double> GetMaxAndMin(std::vector<double> myvector);
      bool findData(const art::Event& evt);
  };

  CosmicFitDisplayNew::CosmicFitDisplayNew(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _mcdiag (conf().mcdiag()),
    _chtag (conf().chtag()),
    _tctag (conf().tctag()),
    _sttag (conf().sttag()),
    doDisplay_ (conf().doDisplay()),
    clickToAdvance_ (conf().clickToAdvance())
    {}

  CosmicFitDisplayNew::~CosmicFitDisplayNew(){}

  void CosmicFitDisplayNew::beginJob() {
    // create diagnostics if requested...
    if ( !doDisplay_ ) return;
    art::ServiceHandle<art::TFileService> tfs;
    directory_ = gDirectory;
    // If needed, create the ROOT interactive environment. See note 1.
    if ( !gApplication ){
    int    tmp_argc(0);
    char** tmp_argv(0);
    application_ = new TApplication( "noapplication", &tmp_argc, tmp_argv );
    }
    // Create a canvas with a guaranteed unique name; the module label is unique within a job.
    TString name  = "canvas_"     + moduleLabel_;
    TString title = "Canvas for " + moduleLabel_;
    int window_size_x(1300);
    int window_size_y(600);
    canvas_ = tfs->make<TCanvas>(name,title,window_size_x,window_size_y);
    canvas_->Divide(2,2);
  }
  void CosmicFitDisplayNew::analyze(const art::Event& event) {
    //Call one or more of the macros from below here, for example:   
    plot2d(event);
  }
     

        
  void CosmicFitDisplayNew::plot2d(const art::Event& event){
    _evt = event.id().event();  
    findData(event);
    std::vector<double> a0, a1, b0, b1, x, y, z;
    //find time clusters:
    unsigned  _ncosmics = _coscol->size();
    unsigned _nch = _chcol->size();
    //loop over tracks:
    for(unsigned int i =0; i < _ncosmics; i++){
      
      CosmicTrackSeed track =(*_coscol)[i];
      if(track._track.converged == false){continue;}
      
      if(isnan(track._track.MinuitParams.A0) == true && isnan(track._track.MinuitParams.A1) == true && isnan(track._track.MinuitParams.B0) == true && isnan(track._track.MinuitParams.B1) == true) continue;

      a0.push_back(track._track.MinuitParams.A0);
      a1.push_back(track._track.MinuitParams.A1);
      b0.push_back(track._track.MinuitParams.B0);
      b1.push_back(track._track.FitParams.B1);
    }

    GeomHandle<Tracker> tracker; 
         
    // Annulus of a cylinder that bounds the tracker/straw info:
    TubsParams envelope(tracker->g4Tracker()->getInnerTrackerEnvelopeParams());
         
    if (doDisplay_) {
              
      std::cout << "Run: " << event.id().run()
      << "  Subrun: " << event.id().subRun()
      << "  Event: " << event.id().event()<<std::endl;
      TLine  major_error_line, minor_error_line, out_line, fit_to_trackxprime, fit_to_trackyprime;
      TArc   arc;
      TPolyMarker poly, out;
      TBox   box;
      TText  text;     
	     
      canvas_->SetTitle("foo title");
      auto pad = canvas_->cd(1);
      pad->Clear();
      canvas_->SetTitle("bar title");
	    double plotLimits(1000.);
      double zlimit{envelope.zHalfLength()};
        
      auto yzplot = pad->DrawFrame(-zlimit,-plotLimits, zlimit, plotLimits);
      yzplot->GetYaxis()->SetTitleOffset(1.25);
      yzplot->SetTitle( "YZ; Z(mm);Y(mm)");
      for(unsigned int i =0; i < _nch; i++){
        ComboHit const& chit =(*_chcol)[i];
			  auto const& p = chit.pos();	
			  auto const& w = chit.wdir();
			  auto const& s = chit.wireRes();
			  //auto const& t = chit.transRes();			
                 
			  /*XYZVec major = (s*w);
			  XYZVec minor = Geom::ZDir().Cross(w) * t;
			  double major2 = (s*w).Mag2();
			  double minor2 = (Geom::ZDir().Cross(w) * t).Mag2();*/
			  
        double z1 = (p.z()+s*w.z());
        double z2 = (p.z()-s*w.z());
        double y1 = (p.y()+s*w.y());
        double y2 = (p.y()-s*w.y());

        major_error_line.DrawLine( z1, y1,z2, y2);
        }
      
              
	    if(a1.size() > 0){
          TF1 *trackline_x = new TF1("line", "[0]+[1]*x", -plotLimits, plotLimits);
          trackline_x->SetParameter(0, b0[0]);
          trackline_x->SetParameter(1, b1[0]);
          trackline_x->SetLineColor(6);
	        //trackline_x->Draw("same");
               
	        pad = canvas_->cd(6);
	        pad->Clear();   

	        auto xyplot = pad->DrawFrame(-plotLimits,-plotLimits,plotLimits,plotLimits);
	        xyplot->GetYaxis()->SetTitleOffset(1.25);
	        xyplot->SetTitle( "XY; X(mm);Y(mm)");
          arc.SetFillStyle(0);
          arc.DrawArc(0.,0., envelope.outerRadius());
          arc.DrawArc(0.,0., envelope.innerRadius());
          for(unsigned int i =0; i < _nch; i++){
            ComboHit const& chit =(*_chcol)[i];
            auto const& p = chit.pos();
            auto const& w = chit.wdir();
            auto const& s = chit.wireRes();	 
            
            double y1 = p.y()+s*w.y();
            double y2 = p.y()-s*w.y();
            double x1 = p.x()+s*w.x();
            double x2 = p.x()-s*w.x();
            major_error_line.DrawLine( x1, y1, x2, y2);
        }
	        
      TF1 *trackline_y = new TF1("line", "[0]+[1]*x", -plotLimits, plotLimits);
      trackline_y->SetParameter(0,b0[0]);
      trackline_y->SetParameter(1, b1[0]);
      trackline_y->SetLineColor(6);
      trackline_y->Draw("same");
      }
      ostringstream title;
      title << "Run: " << event.id().run()
      << "  Subrun: " << event.id().subRun()
      << "  Event: " << event.id().event()<<".root";
    
    text.SetTextAlign(11);
    text.DrawTextNDC( 0., 0.01, title.str().c_str() );

    canvas_->Modified();
    canvas_->Update();
    canvas_->SaveAs(title.str().c_str());
    }
    if ( clickToAdvance_ ){
    	cerr << "Double click in the Canvas " << moduleLabel_ << " to continue:" ;
    	gPad->WaitPrimitive();
  	      } else{
    	char junk;
    	cerr << "Enter any character to continue: ";
    	cin >> junk;
          }
  	        cerr << endl;
  }//End Diag

  void CosmicFitDisplayNew::plot3dXYZ(const art::Event& event){
         
    _evt = event.id().event();  // add event id

    auto comboHits  = event.getValidHandle<ComboHitCollection>( _chtag );
    auto Tracks  = event.getValidHandle<CosmicTrackSeedCollection>( _sttag );
    GeomHandle<Tracker> tracker;

    // Annulus of a cylinder that bounds the tracker
    TubsParams envelope(tracker->g4Tracker()->getInnerTrackerEnvelopeParams());

    std::vector<double> x, y, z, a0, a1, b0, b1;
    std::vector<XYZVec> xprimes, yprimes, zprimes, initial_track_direction;
    x.reserve(comboHits->size());
    y.reserve(comboHits->size());
    z.reserve(comboHits->size());
    a1.reserve(Tracks->size());
    b1.reserve(Tracks->size());
    a0.reserve(Tracks->size());
    b0.reserve(Tracks->size());

    //Get track coordinate system and fit parameters
    for(auto const& track: *Tracks){
      a1.push_back(track._track.MinuitParams.A0);
      a0.push_back(track._track.MinuitParams.A1);
      b1.push_back(track._track.MinuitParams.B0);
      b0.push_back(track._track.MinuitParams.B1);	
    }

 
    for(auto const& chit : *comboHits){
      x.push_back(chit.pos().x());
      y.push_back(chit.pos().y());
      z.push_back(chit.pos().z());
    }

    if (doDisplay_) {
      std::cout << "Run: " << event.id().run()
      << "  Subrun: " << event.id().subRun()
      << "  Event: " << event.id().event()<<std::endl;
      TTUBE tube;
      TText  text;
      TPolyMarker3D poly3D;

      // Draw the frame for the cylinders in plot:
      double plotLimits(1000.);
      double zlimit{envelope.zHalfLength()};
      //Create Pad:
      auto pad = canvas_;
      pad->Clear();
      canvas_->Draw();

      TH3D *xyzplot =new TH3D("XYZ","XYZ",1,-plotLimits,plotLimits,1,-plotLimits,plotLimits,1,-zlimit,zlimit);

      xyzplot->GetYaxis()->SetTitleOffset(1.5);
      xyzplot->GetXaxis()->SetTitleOffset(1.5);
      xyzplot->GetZaxis()->SetTitleOffset(1.5);
      xyzplot->SetTitle( "Visualization in XYZ;X (mm);Y (mm); Z(mm)");
      xyzplot->SetMarkerStyle(2);
      xyzplot->SetMarkerSize(0.65);
      xyzplot->SetMarkerColor(kRed);
      xyzplot->SetStats(0);
      xyzplot->Draw();
      TTUBE *Tube = new TTUBE("","","", envelope.innerRadius(), envelope.outerRadius(), zlimit);
      Tube->SetFillColor(18);
      Tube->SetLineColor(18);
      Tube->Draw();     
      poly3D.SetMarkerStyle(2);
      poly3D.SetMarkerSize(0.65);
      poly3D.SetMarkerColor(kRed);
      int index = 1;
      for ( auto const& chit : *comboHits ){
        TPolyLine3D *errors = new TPolyLine3D;
        errors->SetLineColor(kRed);
        auto const& p = chit.pos();
        auto const& w = chit.wdir();
        auto const& s = chit.wireRes();
        double x1 = p.x()+s*w.x();
        double x2 = p.x()-s*w.x();
        double z1 = p.z()+s*w.z();
        double z2 = p.z()-s*w.z();
        double y1 = p.y()+s*w.y();
        double y2 = p.y()-s*w.y();
        poly3D.SetPoint(index, p.x(), p.y(),p.z() );
        errors->SetPoint(0, x1,y1,z1);
        errors->SetNextPoint(x2,y2,z2);
        index+=1;
        errors->Draw("same");		
      }       
      poly3D.Draw("same");
      if(a1.size() > 0 && b1.size() > 0){
        TPolyLine3D *Track = new TPolyLine3D;
        //TPolyLine3D *InitTrack = new TPolyLine3D;
        Track->SetLineColor(kBlue);
        //InitTrack->SetLineColor(kGreen);  
        double tz1 = -1500;//z[0];
        double tz2 = 1500;//z[z.size()-1];

        double tx1 = a0[0] + a1[0]*tz1;
        double tx2 = a0[0] + a1[0]*tz2;
        double ty1 = b0[0] + b1[0]*tz1;//b1[0]*tz1+b0[0]*tx1+a0[0];
        double ty2 = b0[0] + b1[0]*tz2; //a1[0]*tz2+a0[0]*tx2+b0[0];	
        Track->SetPoint(0, tx1, ty1, tz1);
        Track->SetNextPoint(tx2, ty2, tz2);
        Track->Draw("same");

    }

    canvas_->Update();

    //END 3D
    ostringstream title;


    title << "Run: " << event.id().run()
    << "  Subrun: " << event.id().subRun()
    << "  Event: " << event.id().event()<<".root";


    text.SetTextAlign(11);
    text.DrawTextNDC( 0., 0.01, title.str().c_str() );
    canvas_->Modified();

    canvas_->SaveAs(title.str().c_str());

    if ( clickToAdvance_ ){
    cerr << "Double click in the Canvas " << moduleLabel_ << " to continue:" ;
    gPad->WaitPrimitive();
    } else{
    char junk;
    cerr << "Enter any character to continue: ";
    cin >> junk;
    }
    cerr << endl;
    }

  }//end 3d

      
     
   std::vector<double> CosmicFitDisplayNew::GetMaxAndMin(std::vector<double> myvector){
	std::vector<double> MaxAndMin;
	double first = myvector[0];
    	double smallest = first;
	double biggest = first;
    	for (unsigned i=0; i<  myvector.size(); ++i){ 
	double element = myvector[i];
        if (element< smallest) {
            smallest = element;
         }
    	}
	for (unsigned i=0; i<  myvector.size(); ++i){ 
	double element = myvector[i];
         if (element> biggest) {
            biggest = element;
         }
    	}
    MaxAndMin.push_back(biggest);
    MaxAndMin.push_back(smallest);
    return MaxAndMin;
	}

bool CosmicFitDisplayNew::findData(const art::Event& evt){
	_chcol = 0; 
        _coscol = 0; 
	auto chH = evt.getValidHandle<ComboHitCollection>(_chtag);
	_chcol = chH.product();
	auto stH = evt.getValidHandle<CosmicTrackSeedCollection>(_sttag);
	_coscol =stH.product();
      
	return _chcol != 0 && _coscol !=0 ;
       }   
	
}//End Namespace Mu2e

using mu2e::CosmicFitDisplayNew;
DEFINE_ART_MODULE(CosmicFitDisplayNew);
    

    
