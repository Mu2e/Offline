//Author: S Middleton 
//Date: Jan 2020
//Purpose: To help ana;yse out put of different clustering modules -> run scripts thru: 
//	"valCompare -w validationCalo/ce.html reference.root new_file.root"

//DataProducts:
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
//ROOT
#include "TStyle.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TTree.h"

using namespace std; 

namespace mu2e 
{
  class CaloClusterCompare : public art::EDAnalyzer {
    public:
	struct Config{
		      using Name=fhicl::Name;
		      using Comment=fhicl::Comment;
		      fhicl::Atom<int> diag{Name("diagLevel"), Comment("set to 1 for info"),1};
		      fhicl::Atom<art::InputTag> clustertag{Name("CaloClusterCollection"),Comment("tag for cluster col")};
		      
	       };
	      typedef art::EDAnalyzer::Table<Config> Parameters;

	      explicit CaloClusterCompare(const Parameters& conf);
	      virtual ~CaloClusterCompare();
	      virtual void beginJob() override;
	      virtual void analyze(const art::Event& e) override;
	     
	    private: 
	      
	      	Config _conf;
		int  _diag;
		art::InputTag   _clustertag;

		const CaloClusterCollection* _clustercol;
	      	TH1F* _Energy;
		TH1F* _Time;
		TH1F* _DiskId;
           	TH1F* _EnergyErr;
		TH1F* _TimeErr;
	   	TH1F* _Angle;
	        TH1F* _PosX;  
	    	TH1F* _PosY;
		TH1F* _PosZ; 
  		TH1F* _CryIds;
		TH1F* _ClSize;
		bool findData(const art::Event& evt);
 	};

	CaloClusterCompare::CaloClusterCompare(const Parameters& conf) :
	art::EDAnalyzer(conf),
	_diag (conf().diag()),
	_clustertag (conf().clustertag())
	{}

	CaloClusterCompare::~CaloClusterCompare(){}
	
	void CaloClusterCompare::beginJob() {
      		art::ServiceHandle<art::TFileService> tfs;
		_Energy= tfs->make<TH1F>("Energy Dep [MeV]","EDep" ,120,50, 110);
		_Energy->GetXaxis()->SetTitle("EDep");
		_Time= tfs->make<TH1F>("Time","Time" ,40,500,1700);
		_Time->GetXaxis()->SetTitle("Time");
		_EnergyErr= tfs->make<TH1F>("Energy Dep Err[MeV]","EDep Err" ,20,0, 5);
		_EnergyErr->GetXaxis()->SetTitle("EDepErr");
		_TimeErr= tfs->make<TH1F>("TimeErr","TimeErr" ,50,0,50);
		_TimeErr->GetXaxis()->SetTitle("TimeErr");
		_Angle= tfs->make<TH1F>("Angle","Angle " ,50,-3.1415, 3.1415);
		_Angle->GetXaxis()->SetTitle("Angle");
		_PosX= tfs->make<TH1F>("PosX","Pos X" ,50,-650,650);
		_PosX->GetXaxis()->SetTitle("X Pos");
		_PosY= tfs->make<TH1F>("PosY","PosY " ,50,-650,650);
		_PosY->GetXaxis()->SetTitle("Y Pos");
		_PosZ= tfs->make<TH1F>("PosZ","Pos Z" ,100,1000,1000);
		_PosZ->GetXaxis()->SetTitle("Z Pos");
		_CryIds= tfs->make<TH1F>("Ids","CryIds" ,100,0,2000);
		_CryIds->GetXaxis()->SetTitle("CryIds");
		_ClSize= tfs->make<TH1F>("Clsize","ClSize" ,20,0,20);
		_ClSize->GetXaxis()->SetTitle("ClSize");
		
	}


      	void CaloClusterCompare::analyze(const art::Event& event) {
       
        
		if(!findData(event)) // find data
		throw cet::exception("RECO")<<"No data in event"<< endl; 

		for(size_t ist = 0;ist < _clustercol->size(); ++ist){
		
		CaloCluster sts =(*_clustercol)[ist];
		_Energy->Fill(sts.energyDep());
		_Time->Fill(sts.time());
		_TimeErr->Fill(sts.timeErr());
		_EnergyErr->Fill(sts.energyDepErr());
		_Angle->Fill(sts.angle());
		_PosX->Fill(sts.cog3Vector().x());
		_PosY->Fill(sts.cog3Vector().y());
		_PosZ->Fill(sts.cog3Vector().z());
		_ClSize->Fill(sts.size());
		std::cout<<"========== Event : "<<event.id()<<"========"<<std::endl;
		std::cout<<"Cluster "<<ist<<" size "<<sts.size()<<" Edep "<<sts.energyDep()<<" Time "<<sts.time()<<" Pos "<<sts.cog3Vector().x()<<" "<<sts.cog3Vector().y()<<std::endl;
		for(unsigned i =0 ; i< sts.caloCrystalHitsPtrVector().size();i++){
			art::Ptr< CaloCrystalHit>  cry=sts.caloCrystalHitsPtrVector()[i] ;
           		std::cout<<"crystal ids "<<cry->id()<<" ROid "<<cry->nROId()<<"crystal time "<<cry->time()<<" crystal e dep "<<cry->energyDep()<<" E dep total "<<cry->energyDepTot()<<std::endl;
			
		}
	}
}

	bool CaloClusterCompare::findData(const art::Event& evt){
		_clustercol = 0; 
		
		auto H = evt.getValidHandle<CaloClusterCollection>(_clustertag);
		_clustercol =H.product();
		
		return _clustercol != 0 ;
       }

}

using mu2e::CaloClusterCompare;
DEFINE_ART_MODULE(CaloClusterCompare);



