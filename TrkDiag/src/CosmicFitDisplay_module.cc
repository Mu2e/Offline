
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
#include "DataProducts/inc/threevec.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"

#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/ComboHit.hh"

#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"

#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
// Mu2e Utilities
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
//Mu2e Tracker Geom:
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// Mu2e diagnostics
#include "TrkDiag/inc/ComboHitInfo.hh"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// ROOT incldues

#include "TLegend.h"
#include "TTree.h"
#include "TH2D.h"
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
  class CosmicFitDisplay : public art::EDAnalyzer {
    public:
      explicit CosmicFitDisplay(fhicl::ParameterSet const& pset);
      virtual ~CosmicFitDisplay();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);
    private: 
      //int  _diag;
      bool _mcdiag;
      Int_t _evt; // add event id

      // The module label of this instance of this module.
      std::string moduleLabel_;

      //For Event Displays:
      TApplication* application_;
      TDirectory*   directory_ = nullptr;
      TTree* _cosmic_fit;
      TCanvas*      canvas_ = nullptr;
      TH2D* _display = nullptr;
      TNtuple* _ntTrack = nullptr;
      TNtuple* _ntHit = nullptr;

      // event object Tags

      art::InputTag   _chtag;//combo
      art::InputTag   _tctag;//timeclusters
      art::InputTag   _sttag;//Striaght tracks
      bool doDisplay_;
      bool clickToAdvance_;
      
      
      //Add in details of the tracker segments positions. This will be used for alignment later on:

      Int_t _strawid; // strawid info

      void plot2d(const art::Event& evt);
      void plot3d(const art::Event& evt);
      void improved_event3D(const art::Event& evt);
      std::vector<double> GetMaxAndMin(std::vector<double> myvector);
    };
    CosmicFitDisplay::CosmicFitDisplay(fhicl::ParameterSet const& pset) :
	art::EDAnalyzer(pset),
	
	_mcdiag		(pset.get<bool>("MCdiag",true)),
	_chtag		(pset.get<art::InputTag>("ComboHitCollection")),
	_tctag		(pset.get<art::InputTag>("TimeClusterCollection")),
	_sttag		(pset.get<art::InputTag>("CosmicTrackSeedCollection")),
	doDisplay_(pset.get<bool>("doDisplay",false)),
        clickToAdvance_(pset.get<bool>("clickToAdvance",false))
        {}

    CosmicFitDisplay::~CosmicFitDisplay(){}

    void CosmicFitDisplay::beginJob() {
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



      void CosmicFitDisplay::analyze(const art::Event& event) {
        plot2d( event);
        
      }//End Analyze 


  
      void CosmicFitDisplay::plot2d(const art::Event& event){
        
        _evt = event.id().event();  
       
        auto comboHits  = event.getValidHandle<ComboHitCollection>( _chtag );
        auto Tracks  = event.getValidHandle<CosmicTrackSeedCollection>( _sttag );

        std::vector<double> x, y, z, a0, a1, b0,b1;
        x.reserve(comboHits->size());
        y.reserve(comboHits->size());
        z.reserve(comboHits->size());
	
        a0.reserve(Tracks->size());
        a1.reserve(Tracks->size());
        b0.reserve(Tracks->size());
        b1.reserve(Tracks->size());
        
        // loop over combo hits
        for(auto const& chit : *comboHits){
        x.push_back(chit.pos().x());
        y.push_back(chit.pos().y());
        z.push_back(chit.pos().z());
        }

        //loop over tracks:
        for(auto const& track: *Tracks){
        a0.push_back(track._track.get_track_parameters()[0]);
        a1.push_back(track._track.get_track_parameters()[1]);
        b0.push_back(track._track.get_track_parameters()[2]);
        b1.push_back(track._track.get_track_parameters()[3]);
        }
        
        GeomHandle<Tracker> tracker;
        //Straw& straw_details = tracker->getStraw(id);
        // Annulus of a cylinder that bounds the tracker/straw info:
        TubsParams envelope(tracker->getInnerTrackerEnvelopeParams());
	//double const straw_radius = straw_details.at(0).getRadius();; 
        

        if (doDisplay_) {
              std::cout << "Run: " << event.id().run()
           << "  Subrun: " << event.id().subRun()
           << "  Event: " << event.id().event()<<std::endl;
              TLine  line, fit_to_track;
	      TArc   arc;
              TPolyMarker poly;
	      TBox   box;
	      TText  text;
	      
	      arc.SetFillStyle(0);
	      //straw.SetFillStyle(0);
	      //poly.SetMarkerStyle(2);
	      //poly.SetMarkerSize(straw_radius);
	      //poly.SetMarkerColor(kBlue);

	      canvas_->SetTitle("foo title");
	      auto pad = canvas_->cd(1);
	      pad->Clear();
	      canvas_->SetTitle("bar title");

	      // Draw the frame for the y vs x plot.
	      double plotLimits(1000.);
	      auto xyplot = pad->DrawFrame(-plotLimits,-plotLimits,plotLimits,plotLimits);

	      // Draw the inner and outer arc of the tracker.
	      arc.SetLineColor(kBlack);
	      arc.DrawArc(0.,0., envelope.outerRadius());
	      arc.DrawArc(0.,0., envelope.innerRadius());

	      xyplot->GetYaxis()->SetTitleOffset(1.25);
	      xyplot->SetTitle( "y vs x;(mm);(mm)");

	      line.SetLineColor(kRed);
	      fit_to_track.SetLineColor(kBlue);
	      poly.SetMarkerColor(kRed);
              
              for ( auto const& chit : *comboHits ){
                        
			auto const& p = chit.pos();
			auto const& w = chit.wdir();
			auto const& s = chit.wireRes();
			double x0{p.x()};
			double y0{p.z()};
                        //double z0{p.z()};
			poly.DrawPolyMarker( 1, &x0, &y0 );
		        
			double x1 = p.x()+s*w.x();
			double x2 = p.x()-s*w.x();
			double y1 = p.y()+s*w.y();
			double y2 = p.y()-s*w.y();
			line.DrawLine( x1, y1, x2, y2);
                  
      	      }
 	   
              
	      if(a1.size() > 0){
		
		double tx1 = z[0];
		double tx2 = z[z.size()-1];
		double ty1 = a1[0]*tx1+a0[0];
		double ty2 =a1[0]*tx2+a0[0];
		fit_to_track.DrawLine( tx1, ty1, tx2, ty2);
		 

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
      
      }//End 2d

      void CosmicFitDisplay::plot3d(const art::Event& event){
         // find data in event
        //findData(event);
        _evt = event.id().event();  // add event id

        //get combo hits
        auto comboHits  = event.getValidHandle<ComboHitCollection>( _chtag );
        auto Tracks  = event.getValidHandle<CosmicTrackSeedCollection>( _sttag );
        GeomHandle<Tracker> tracker;

        // Annulus of a cylinder that bounds the tracker
        TubsParams envelope(tracker->getInnerTrackerEnvelopeParams());


	// Create arrays for x,y,z:
	std::vector<double> x, y, z, a0, a1, b0, b1;
        x.reserve(comboHits->size());
        y.reserve(comboHits->size());
        z.reserve(comboHits->size());

        a1.reserve(Tracks->size());
	b1.reserve(Tracks->size());
        a0.reserve(Tracks->size());
        b0.reserve(Tracks->size());
       
        // loop over combo hits
        for(auto const& chit : *comboHits){
        x.push_back(chit.pos().x());
        y.push_back(chit.pos().y());
        z.push_back(chit.pos().z());
        }
        //loop over tracks:
        for(auto const& track: *Tracks){
        a1.push_back(track._track.get_track_parameters()[1]);
        a0.push_back(track._track.get_track_parameters()[0]);
        b1.push_back(track._track.get_track_parameters()[3]);
        b0.push_back(track._track.get_track_parameters()[2]);
	//mz.push_back(track._track.get_m_1());
        }
        
        //If Diag:
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
	      xyzplot->SetTitle( "Visualization in XYZ;x (mm);y (mm); z(mm)");
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
			double x0{p.x()};
		        double y0{p.y()};
                        double z0{p.z()};
                        auto const& w = chit.wdir();
		        auto const& s = chit.wireRes();
			double x1 = p.x()+s*w.x();
			double x2 = p.x()-s*w.x();
                        double z1 = p.z()+s*w.z();
        		double z2 = p.z()-s*w.z();
        		double y1 = p.y()+s*w.y();
        		double y2 = p.y()-s*w.y();
			poly3D.SetPoint(index, x0, y0, z0 );
			errors->SetPoint(0, x1,y1,z1);
			errors->SetNextPoint(x2,y2,z2);
		        index+=1;
                        errors->Draw("same");
			
      	      }
              
              poly3D.Draw("same");

	      if(a1.size() > 0 && b1.size() > 0){
	      	TPolyLine3D *Track = new TPolyLine3D;
                Track->SetLineColor(kBlue);
                double tx1 = x[0];
	        double tx2 = x[x.size()-1];
		double tz1 = z[0];
	        double tz2 = z[z.size()-1];
		double ty1 = b1[0]*tz1+b0[0]*tx1+a0[0];
		double ty2 =a1[0]*tz2+a0[0]*tx2+b0[0];
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

           }//end display
      }//end 3d

      
      void CosmicFitDisplay::improved_event3D(const art::Event& event){
          
          TGeoManager *geom = new TGeoManager("tracker_geom", "tracker_geom");
          //Materials/Media
          TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
          TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
          //Set Top Geom:
          TGeoVolume *top = geom->MakeBox("3DGeometry", Vacuum, 2000,2000,2000); 
           geom->SetTopVolume(top);
   
           TGeoVolume *rootbox = geom->MakeBox("ROOT", Vacuum, 110., 50., 5.);
           rootbox->SetVisibility(kFALSE);
          //Rotate about z
          TGeoRotation *rot = new TGeoRotation("rot",0,90,0.);
          rot->RegisterYourself();

          GeomHandle<Tracker> tracker;
 
          // Annulus of a cylinder that bounds the tracker
          TubsParams envelope(tracker->getInnerTrackerEnvelopeParams());
          
	  // Draw the frame for the cylinders in plot:
          double zlimit{envelope.zHalfLength()};
          double R_min{envelope.innerRadius()};
          double R_max{envelope.outerRadius()};
         
          
          TGeoVolume *tr = geom->MakeTube("tr", Vacuum, R_min, R_max, zlimit);
	   tr->SetLineColor(kViolet);

_evt = event.id().event();  // add event id

            //get combo hits
           auto comboHits  = event.getValidHandle<ComboHitCollection>( _chtag );
           auto Tracks  = event.getValidHandle<CosmicTrackSeedCollection>( _sttag );
        


	   // Create arrays for x,y,z:
	   std::vector<double> x, y, z,a0,a1,b0,b1;
           x.reserve(comboHits->size());
           y.reserve(comboHits->size());
           z.reserve(comboHits->size());

          a1.reserve(Tracks->size());
	  b1.reserve(Tracks->size());
          a0.reserve(Tracks->size());
          b0.reserve(Tracks->size());
          // loop over combo hits
          for(auto const& chit : *comboHits){
             x.push_back(chit.pos().x());
             y.push_back(chit.pos().y());
             z.push_back(chit.pos().z());
          }
          //loop over tracks:
          for(auto const& track: *Tracks){
             a1.push_back(track._track.get_track_parameters()[1]);
             a0.push_back(track._track.get_track_parameters()[0]);
             
	     b1.push_back(track._track.get_track_parameters()[3]);
             b1.push_back(track._track.get_track_parameters()[2]);
          }
        
        //If Diag:
        if (doDisplay_) {
              std::cout << "Run: " << event.id().run()
           << "  Subrun: " << event.id().subRun()
           << "  Event: " << event.id().event()<<std::endl;


           //Add Node
	   top->AddNode(tr,1, rot);
          

          //Close the geometry
          //geom->CloseGeometry();
          geom->SetVisLevel(4);
          top->Draw("ogle");
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
   }

   std::vector<double> CosmicFitDisplay::GetMaxAndMin(std::vector<double> myvector){
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

   
	
}//End Namespace Mu2e

using mu2e::CosmicFitDisplay;
DEFINE_ART_MODULE(CosmicFitDisplay);
    

    


