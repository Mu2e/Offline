//Author: S Middleton
//Date: Nov 2019
//Purpose: For the purpose of fast crystal hit from recodigi

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"

#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <string>
#include <cmath>

namespace mu2e {
	class FastCaloCrystalHitFromHit : public art::EDProducer {
  		public:

			explicit FastCaloCrystalHitFromHit(fhicl::ParameterSet const& pset) :
			art::EDProducer{pset},
			caloDigisToken_{consumes<CaloRecoDigiCollection>(pset.get<std::string>("caloDigisModuleLabel"))},
			time4Merge_          (pset.get<double>     ("time4Merge")),
			diagLevel_           (pset.get<int>        ("diagLevel",0))
			{
				produces<CaloCrystalHitCollection>();
			}

	    		void produce(art::Event& e) override;

  		private:

			typedef art::Ptr<CaloRecoDigi> CaloRecoDigiPtr;

			art::ProductToken<CaloRecoDigiCollection> const caloDigisToken_;
			double      time4Merge_;
			int         diagLevel_;
			std::vector<std::vector<const CaloRecoDigi*>> hitMap_;  
			void MakeCaloCrystalHits(CaloCrystalHitCollection& CaloCrystalHits, const art::ValidHandle<CaloRecoDigiCollection>& recoCaloDigisHandle);
			void FillBuffer(int crystalId, int nRoid, double time, double timeErr, double eDep, double eDepErr,
                    std::vector<CaloRecoDigiPtr>& buffer, CaloCrystalHitCollection& CaloCrystalHits);
	};

	void FastCaloCrystalHitFromHit::produce(art::Event& event)
  	{
		auto const& recoCaloDigisHandle = event.getValidHandle(caloDigisToken_);
		auto CaloCrystalHits = std::make_unique<CaloCrystalHitCollection>();
		MakeCaloCrystalHits(*CaloCrystalHits, recoCaloDigisHandle);
		//store the crystal hit collection:
		event.put(std::move(CaloCrystalHits));
 	 }

	void FastCaloCrystalHitFromHit::MakeCaloCrystalHits(CaloCrystalHitCollection& CaloCrystalHits, art::ValidHandle<CaloRecoDigiCollection> const& recoCaloDigisHandle)
  	{
		Calorimeter const &cal = *(GeomHandle<Calorimeter>());
		auto const& recoCaloDigis = *recoCaloDigisHandle;
		CaloRecoDigi const* base = &recoCaloDigis.front(); // TODO What if recoCaloDigis is empty?
		if (cal.nRO() > int(hitMap_.size())) continue;
		
		for (int i = hitMap_.size(); i<= cal.nRO(); ++i) hitMap_.push_back(std::vector<const CaloRecoDigi*>());
		for (size_t i=0; i<hitMap_.size(); ++i) hitMap_[i].clear();

		for (unsigned int i=0; i< recoCaloDigis.size(); ++i)
		{
			int crystalId = cal.caloInfo().crystalByRO(recoCaloDigis[i].ROid());
			hitMap_[crystalId].push_back(&recoCaloDigis[i]);
		}


		for (unsigned int crystalId=0;crystalId<hitMap_.size();++crystalId)
		{
			std::vector<const CaloRecoDigi*> &hits = hitMap_[crystalId];
        
			//check if empty:
			if (hits.empty()) continue;

			//sort hits in terms of time:
			std::sort(hits.begin(),hits.end(),[](const auto a, const auto b){return a->time() < b->time();});

			//find hit:
			auto startHit = hits.begin();
			auto endHit   = hits.begin();

			//create a buffer 
			std::vector<CaloRecoDigiPtr> buffer;
			double timeW(0);
			double eDepTot(0),eDepTotErr(0);
			int nRoid(0);

			//loop through hits:
			while (endHit != hits.end())
          		{
				//time:
				double deltaTime = (*endHit)->time()-(*startHit)->time();

				//if > than set merge time:
				if (deltaTime > time4Merge_) 
				{
					double time = timeW/nRoid;
					double timeErr = 0;
					//fill that buffer: //TODO -->peakpos, flags etc.
					FillBuffer(crystalId, nRoid, time, timeErr, eDepTot/nRoid, eDepTotErr/nRoid, buffer, CaloCrystalHits);
					//clear:
					buffer.clear();
					timeW      = 0.0;
					eDepTot    = 0.0;
					eDepTotErr = 0.0;
					nRoid      = 0;
					startHit   = endHit;
              			}
            			else //if not then:
              			{
					//add up times
					timeW      += (*endHit)->time();
					//add up energies
					eDepTot    += (*endHit)->energyDep();
					eDepTotErr += (*endHit)->energyDepErr() * (*endHit)->energyDepErr();
					//increment no. ROs
					++nRoid;
					//get index:
					size_t index = *endHit - base;
					//add to end of hit buffer:
					buffer.push_back(art::Ptr<CaloRecoDigi>(recoCaloDigisHandle, index));
					//increment:
					++endHit;
              			}

          		}
			//get time and error
			double time = timeW/nRoid;
			double timeErr = 0;
			//Finally fill buffer:
			FillBuffer(crystalId, nRoid, time, timeErr, eDepTot/nRoid, eDepTotErr/nRoid, buffer, CaloCrystalHits);
		}


		if ( diagLevel_ > 0 )
		{
			printf("[FastCaloCrystalHitFromHit::produce] produced RecoCrystalHits ");
			printf(": CaloCrystalHits.size()  = %i \n", int(CaloCrystalHits.size()));
		}

	}

	void FastCaloCrystalHitFromHit::FillBuffer(int const crystalId,
                                         int const nRoid,
                                         double const time,
                                         double const timeErr,
                                         double const eDep,
                                         double const eDepErr,
                                         std::vector<CaloRecoDigiPtr>& buffer,
                                         CaloCrystalHitCollection& CaloCrystalHits)
  	{
 
	    CaloCrystalHits.emplace_back(CaloCrystalHit(crystalId, nRoid, time, timeErr, eDep, eDepErr, buffer));

	    if (diagLevel_ > 1)
	      {
			std::cout<<"[FastCaloCrystalHitFromHit] created hit in crystal id="<<crystalId<<"\t with time="<<time<<"\t eDep="<<eDep<<"\t  from "<<nRoid<<" RO"<<std::endl;

		
	      }
  	}	


}

DEFINE_ART_MODULE(mu2e::FastCaloCrystalHitFromHit);
