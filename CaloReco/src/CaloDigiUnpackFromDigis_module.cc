//
// Recall, the caloFigiColl structure is the following:
// nTotWords - nWords_roID - roiID - nWord_roID_ihit - time_roID_ihit - Data_roID_ihit - ...
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloDigiPacked.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"



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


       void makeCaloRecoDigis(const CaloDigiPacked& caloDigis, CaloDigiCollection& caloDigiColl);       

  };



  //-------------------------------------------------------
  void CaloDigiUnpackFromDigis::produce(art::Event& event) 
  {

       if (diagLevel_ > 0 ) printf("[CaloDigiUnpackFromDigis::produce] event %i\n", event.id().event());

       //Get handles to calorimeter RO (aka APD) collection
       art::Handle<CaloDigiPacked> caloDigisHandle;
       event.getByLabel(caloDigiModuleLabel_, caloDigisHandle);
       const CaloDigiPacked& caloDigis(*caloDigisHandle);


       std::unique_ptr<CaloDigiCollection> caloDigiColl(new CaloDigiCollection);
       makeCaloRecoDigis(caloDigis, *caloDigiColl);       
       event.put(std::move(caloDigiColl));

       return;
  }

    
  //--------------------------------------------------------------------------------------
  void CaloDigiUnpackFromDigis::makeCaloRecoDigis(const CaloDigiPacked& caloDigis, CaloDigiCollection& caloDigiColl)
                                                   
  {
     
      std::vector<int> caloFromDigi = caloDigis.output();
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
