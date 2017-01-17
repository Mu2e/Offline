//
// Recall, the caloFigiColl structure is the following:
// nTotWords - nWords_roID - roiID - nWord_roID_ihit - time_roID_ihit - Data_roID_ihit - ...
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "RecoDataProducts/inc/CaloDigiPackedCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"

#include <iostream>
#include <string>


namespace mu2e {



  class CaloDigiUnpackFromDigis : public art::EDProducer {

    public:

        explicit CaloDigiUnpackFromDigis(fhicl::ParameterSet const& pset) :
          caloDigiModuleLabel_ (pset.get<std::string>("caloDigiModuleLabel")),
          minDigiHitLength_    (pset.get<int>        ("minDigiHitLength")),
          diagLevel_           (pset.get<int>        ("diagLevel",0))
        {

            produces<CaloDigiCollection>();
        }

        virtual ~CaloDigiUnpackFromDigis() {}
        void produce(art::Event& e);



    private:

       std::string  caloDigiModuleLabel_; 
       int          minDigiHitLength_;
       int          diagLevel_;      
       double       edepTot_;


       void makeCaloRecoDigis(const CaloDigiPackedCollection& caloDigisPacked, CaloDigiCollection& caloDigiColl);       
       void fillDigis(const CaloDigiPacked& caloDigiPacked, CaloDigiCollection& caloDigiColl);       

  };



  //-------------------------------------------------------
  void CaloDigiUnpackFromDigis::produce(art::Event& event) 
  {

       if (diagLevel_ > 0 ) printf("[CaloDigiUnpackFromDigis::produce] start event %i\n", event.id().event());

       //Get handles to calorimeter RO (aka APD) collection
       art::Handle<CaloDigiPackedCollection> caloDigisHandle;
       event.getByLabel(caloDigiModuleLabel_, caloDigisHandle);
       const CaloDigiPackedCollection& caloDigisPacked(*caloDigisHandle);

       std::unique_ptr<CaloDigiCollection> caloDigiColl(new CaloDigiCollection);              
       makeCaloRecoDigis(caloDigisPacked, *caloDigiColl);                     
       event.put(std::move(caloDigiColl));
    
       if ( diagLevel_ > 0 ) std::cout << "[CaloDigiUnpackFromDigis::produce] end" << std::endl;

       return;
  }

    
  //--------------------------------------------------------------------------------------
  void CaloDigiUnpackFromDigis::makeCaloRecoDigis(const CaloDigiPackedCollection& caloDigisPacked, CaloDigiCollection& caloDigiColl)
                                                   
  {
     for (const auto& caloDigiPacked : caloDigisPacked) fillDigis(caloDigiPacked,caloDigiColl);
  }          


  //--------------------------------------------------------------------------------------
  void CaloDigiUnpackFromDigis::fillDigis(const CaloDigiPacked& caloDigiPacked, CaloDigiCollection& caloDigiColl)
                                                   
  {     
      std::vector<int> caloFromDigi = caloDigiPacked.output();
      unsigned int caloFromDigiSize = caloFromDigi.size();

      unsigned int index(1); //starts form one, the element 0 is the size of the CaloDigi array
      while ( index < caloFromDigiSize )
      {        
           int    digitizedHitLength = caloFromDigi.at(index);
           int    roId               = caloFromDigi.at(index+1);

           if (digitizedHitLength > minDigiHitLength_) 
           {               
               int wfSize   = digitizedHitLength - 2;
               int wfOffset = index+2;
               
               std::vector<int> waveform(&caloFromDigi[wfOffset],&caloFromDigi[wfOffset+wfSize]);
              
               unsigned int idx(0);
               while(idx < waveform.size())
               {        
                   int    size    = waveform.at(idx);
                   double t0      = waveform.at(idx+1);
                   int    wfsize  = size-2;
                   int    wfindex = idx+2;
                                      
                   std::vector<int> sub(&waveform[wfindex],&waveform[wfindex+wfsize]);
                   caloDigiColl.emplace_back( CaloDigi(roId,t0,sub) );         
                   idx += size;
               }
           }
           
           index += digitizedHitLength;           
      }
      
      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiUnpackFromDigis] caloDigiSize="<<caloDigiColl.size()<<std::endl;
  }          


}

using mu2e::CaloDigiUnpackFromDigis;
DEFINE_ART_MODULE(CaloDigiUnpackFromDigis);
