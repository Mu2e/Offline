// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <list>
#include <deque>
//Mu2e Geom:
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"

//Mu2e Data Prods:
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
// Mu2e Utilities
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"

// Mu2e diagnostics
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// BaBar Kalman filter includes
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/BaBar/BaBar.hh"

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
  class CaloEventDisplays : public art::EDAnalyzer {
    public:
	    struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;

        fhicl::Atom<art::InputTag> calocrysTag{Name("CaloCrystalHitCollection"),Comment("cal reco crystal hit info")};

        fhicl::Atom<art::InputTag> caloclusterTag{Name("CaloClusterCollection"),Comment("cal reco cluster info")};
        
        fhicl::Atom<bool> doDisplay{Name("doDisplay"),Comment("use display"), false};
        fhicl::Atom<bool> clickToAdvance{Name("clickToAdvance"),Comment("next event"), false};
	     };
	    typedef art::EDAnalyzer::Table<Config> Parameters;

          	explicit CaloEventDisplays(const Parameters& conf);
	    
          	virtual ~CaloEventDisplays();
          	virtual void beginJob();
          	virtual void analyze(const art::Event& e) override;
    private: 
	    Config _conf;
    	bool _mcdiag;
    	Int_t _evt; 

	    // The module label of this instance of this module.
	    std::string moduleLabel_;
	    
	    //For Event Displays:
	    TApplication* application_;
	    TDirectory*   directory_ = nullptr;
	    TCanvas*      canvas_ = nullptr;
	    TH2D* _display = nullptr;
	    TNtuple* _ntTrack = nullptr;
	    TNtuple* _ntHit = nullptr;

	    
      art::InputTag _calocrysTag;
      art::InputTag _caloclusterTag;
      art::InputTag _kalrepTag;
      const CaloCrystalHitCollection*  _calcryhitcol;
      const CaloClusterCollection* _calclustercol;


      bool doDisplay_;
      bool clickToAdvance_;
      void plot2d(const art::Event& evt);

      bool findData(const art::Event& evt);
    };

    CaloEventDisplays::CaloEventDisplays(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _calocrysTag(conf().calocrysTag()),
    _caloclusterTag(conf().caloclusterTag()),
    doDisplay_ (conf().doDisplay()),
    clickToAdvance_ (conf().clickToAdvance())
    {}
   
    CaloEventDisplays::~CaloEventDisplays(){}

    void CaloEventDisplays::beginJob() {
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
      //canvas_->Divide(1,2);
   
  }

  void CaloEventDisplays::analyze(const art::Event& event) {

	  plot2d(event);
  }

  void CaloEventDisplays::plot2d(const art::Event& event){
  _evt = event.id().event();  
  findData(event);

  std::vector<double>  clusterxs;

  //unsigned  _ncluster = _calclustercol->size();
  unsigned _ncrystalhits = _calcryhitcol->size();

  art::ServiceHandle<GeometryService> geom;
  if( ! geom->hasElement<Calorimeter>() ) return;
  Calorimeter const & cal = *(GeomHandle<Calorimeter>());
	
	if (doDisplay_) {
      
    std::cout << "Run: " << event.id().run()
    << "  Subrun: " << event.id().subRun()
    << "  Event: " << event.id().event()<<std::endl;


    TText  text;     
    TArc arcOut, arcIn;
    TBox box;
    canvas_->SetTitle("foo title");
    auto pad = canvas_->cd();
    pad->Clear();
    canvas_->SetTitle("bar title");

    auto xyplot = pad->DrawFrame(-1000,-1000, 1000,1000);
    xyplot->GetYaxis()->SetTitleOffset(1.25);
    xyplot->SetTitle( "View of Calo Disk 1 in YZ Plane; Z(mm);Y(mm)");

    art::ServiceHandle<mu2e::GeometryService>   geom;

    mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
    mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;

    float _clusterEdep = 0;
    for (unsigned int tclu=0; tclu<_calclustercol->size();++tclu){
    CaloCluster const& cluster = (*_calclustercol)[tclu];
    _clusterEdep      = cluster.energyDep();
    }

    Disk const & disk =  cal.disk(1);
    double outerR = disk.outerRadius();
    double innerR= disk.innerRadius();

    arcOut.SetFillColor(kGray);
    arcIn.SetFillColor(kWhite);
    arcIn.SetLineColor(kGray+1);
    arcOut.SetLineColor(kGray+1);
    arcOut.DrawArc(0.,0., outerR);
    arcIn.DrawArc(0.,0., innerR);

    typedef std::list<CaloCrystalHit const*> CaloCrystalList;
    CaloCrystalList firstlist;
    for(size_t i =0; i < _ncrystalhits; i++){ 
      CaloCrystalHit const& hit =(*_calcryhitcol)[i];
      firstlist.push_back(&hit);
     }

    firstlist.sort([] (CaloCrystalHit const* lhs, CaloCrystalHit const* rhs) {return lhs->energyDep() > rhs->energyDep();} );

    unsigned int i =0;
    for(int i=0;i<674;i++){
	    Crystal const &crystal = cal.crystal(i);
      double crystalXLen = crystal.size().x();
      double crystalYLen = crystal.size().y();
      CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(1,crystal.position());
      box.SetLineColor(kGray+1);
      box.DrawBox(crystalPos.x()-crystalXLen/2, crystalPos.y()-crystalYLen/2,crystalPos.x()+crystalXLen/2, crystalPos.y()+crystalYLen/2);
		
	  }
    for(auto const& hit : firstlist){
      //CaloCrystalHit const& hit =(*_calcryhitcol).at(i);//_calcryhitcol
      int diskId     = cal.crystal(hit->id()).diskId();

      CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit->id()).position());
      int crystalID = hit->id();
      Crystal const &crystal = cal.crystal(crystalID);
      double crystalXLen = crystal.size().x();
      double crystalYLen = crystal.size().y();
      if(i==0) box.SetFillColor(kRed);
      if(i==1 or i==2) box.SetFillColor(kOrange);
      if(i==3 or i==4) box.SetFillColor(kYellow);
      if(i>5 and i < 7) box.SetFillColor(kGreen);
      if(i>7) box.SetFillColor(kCyan);
      box.DrawBox(crystalPos.x()-crystalXLen/2, crystalPos.y()-crystalYLen/2,crystalPos.x()+crystalXLen/2, crystalPos.y()+crystalYLen/2);

      if(hit->energyDep()>0 and _clusterEdep>0){
					   
        TLatex latex;
        stringstream crys;
        crys<<round(hit->energyDep());
        const char* str_crys = crys.str().c_str();
        latex.SetTextSize(0.02);
        latex.DrawLatex(crystalPos.x()-crystalXLen/2, crystalPos.y()-crystalYLen/2,str_crys);
        cout<<i<<" "<<round(hit->energyDep())<<endl;
        }
        i++;
  	  }
	            
       
      ostringstream title;
      title << "Run: " << event.id().run()
      << "  Subrun: " << event.id().subRun()
      << "  Event: " << event.id().event()<<".root";

      text.SetTextAlign(11);
      text.DrawTextNDC( 0., 0.01, title.str().c_str());

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
    }
	}//display
    
 

  bool CaloEventDisplays::findData(const art::Event& evt){

	
	  _calcryhitcol =0;
	  _calclustercol=0;
	
    auto cryhit = evt.getValidHandle<CaloCrystalHitCollection>(_calocrysTag);
    _calcryhitcol =cryhit.product();
    auto cluster= evt.getValidHandle<CaloClusterCollection>(_caloclusterTag);
    _calclustercol =cluster.product();
        
	  return _calcryhitcol!=0 && _calclustercol !=0 ;
  }


}  // end namespace mu2e 

using mu2e::CaloEventDisplays;
DEFINE_ART_MODULE(CaloEventDisplays);

