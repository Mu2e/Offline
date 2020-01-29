//Author: S Middleton 
//Date: Jan 2020
//Purpose: To help ana;yse out put of different clustering modules

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
		_Energy= tfs->make<TH1F>("Energy Dep [MeV]","EDep " ,50,0, 110);
		_Energy->GetXaxis()->SetTitle("EDep");
		_Energy->SetStats();
	}


      	void CaloClusterCompare::analyze(const art::Event& event) {
       
        
		if(!findData(event)) // find data
		throw cet::exception("RECO")<<"No Time Clusters in event"<< endl; 

		for(size_t ist = 0;ist < _clustercol->size(); ++ist){
		CaloCluster sts =(*_clustercol)[ist];
		_Energy->Fill(sts.energyDep());
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



