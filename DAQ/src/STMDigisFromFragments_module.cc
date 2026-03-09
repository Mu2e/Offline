// =====================================================================
//
// STMDigisFromFragments: create all types of STMDigis from STMFragments
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMMWDDigi.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>

#include <string>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <list>
#include <numeric>
#include <random>
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TLine.h>
#include <TGraph.h>
#include <stdbool.h> //temporary

namespace art
{
  class STMDigisFromFragments;
}

using art::STMDigisFromFragments;

// ======================================================================

class art::STMDigisFromFragments : public EDProducer
{
  public:
  struct Config
  {
    fhicl::Atom<art::InputTag> stmTag {fhicl::Name("stmTag"), fhicl::Comment("Input module")};
    // TODO: add fhicl parameters to that we can choose which types of fragments we read out
    // e.g.     fhicl::Atom<bool> processRaw {fhicl::Name("processRaw"), fhicl::Comment("Process Raw STMFragments")};
  };

  // --- C'tor/d'tor:
  explicit STMDigisFromFragments(const art::EDProducer::Table<Config>& config);

  // --- Production:
  virtual void produce(Event&);

  private:

  art::InputTag _stmFragmentsTag;

}; // STMDigisFromFragments

// ======================================================================


STMDigisFromFragments::STMDigisFromFragments(const art::EDProducer::Table<Config>& config) :
  art::EDProducer{config}
  ,_stmFragmentsTag(config().stmTag())

{
  // Set the size of the vector
  produces<mu2e::STMWaveformDigiCollection>("raw");
  produces<mu2e::STMWaveformDigiCollection>("zs");
  produces<mu2e::STMMWDDigiCollection>("mwd"); // TODO: we should create an STMMWDDigi collection instead of STMWaveformDigis for this
}

// ----------------------------------------------------------------------


void STMDigisFromFragments::produce(Event& event)
{
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> zs_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMMWDDigiCollection> mwd_digis(new mu2e::STMMWDDigiCollection);
  // std::unique_ptr<mu2e::STMWaveformDigiCollection> mwd_waveform_digis(new mu2e::STMWaveformDigiCollection);//Original

  art::Handle<artdaq::Fragments> STMFragmentsH;
  event.getByLabel(_stmFragmentsTag, STMFragmentsH);
  const auto STMContainerFragments = STMFragmentsH.product();


    //Variables to keep track of all zero element arrays& empty arrays
  
    size_t zeroRaw_arrays = 0;
    size_t zeroZS_arrays = 0;
    size_t emptyRaw_arrays = 0;
    size_t emptyZS_arrays = 0;
    
  for (const auto& frag : *STMContainerFragments) {
    // const auto& stm_frag = static_cast<mu2e::STMFragment>(frag); //What was here originally to read other file
    artdaq::ContainerFragment contf(frag);

    
    for (size_t ii = 0; ii < contf.block_count(); ++ii){
      const artdaq::Fragment art_frag = *contf[ii];
      const auto& stm_frag = static_cast<mu2e::STMFragment>(art_frag);
      mu2e::STMWaveformDigi stm_waveform;

      // std::cout<<"ii = "<< ii<<" isRaw = "<< stm_frag.isRaw()<<" isZS = "<<stm_frag.isZS()<<" isMWD = "<<stm_frag.isMWD()<<std::endl;//Check for what we are reading per index
      
       if (stm_frag.isRaw()) {
	 
	 //A check to see what inside these arrays
	 auto payloadarray= stm_frag.payloadBegin();
	 auto N = stm_frag.payloadWords();	 
	 bool allZeros = true; //Assume we have an array with all zeros
	 bool first20hasZeros = false; //Checking header
	 
	 if(N ==0){ //If payloadwords = 0 then the array is empty
	   ++emptyRaw_arrays;
	   std::cout<< "Found empty Raw array at ii = "<<ii<< " \n";// extra check ->Raw has no empty arrays it seems
	 }else{
	   for(size_t k = 0; k< N; ++k){
	     if(payloadarray[k] != 0){
	       allZeros = false; //Break if you find an element that isn't zero
	       break;
	     }
	   }
	   //Check to see if first 20 elements have any zero
	   for(size_t z20 =0; z20 < std::min<size_t>(N,20); ++z20){
	     if(payloadarray[z20] == 0){
	       first20hasZeros = true;//Break if you find any element with a zero, will disrupt header
	       break;
	     }
	   }
	   
	   if (allZeros == true || first20hasZeros == true){
	     ++zeroRaw_arrays; //Adds to the zero array counter
	     std::cout<<"PayloadWords for zeroRaw array: "<<N<< std::endl; //extra check
	     for(size_t zz = 0; zz < std::min<size_t>(N,10);++zz){
	       std::cout<<payloadarray[zz]<<" "; //Prints out first 10 elemenmts to check we have all zeros 
	     }
	     std::cout<<" Raw ii = " << ii<< "\n";
	       
	   } else {
	     for(size_t kk = 0; kk < std::min<size_t>(N,20); ++kk){
	       std::cout<<payloadarray[kk]<<" "; // prints out the first 20 elements
	     }
	     std::cout<<" Raw ii = "<< ii <<"\n";
	     stm_waveform.set_data(stm_frag.payloadWords(), stm_frag.payloadBegin());
	     raw_waveform_digis->emplace_back(stm_waveform);
	   }
	 }
       
	 //end of check

	 //Originally outside the check
	 //stm_waveform.set_data(stm_frag.payloadWords(), stm_frag.payloadBegin());
	 //raw_waveform_digis->emplace_back(stm_waveform);
	 
       }
       else if (stm_frag.isZS()) {
         //A check to see what inside these arrays                                                                                                                                                                             
	 auto payloadarray= stm_frag.payloadBegin();
         auto N = stm_frag.payloadWords();
         bool allZeros = true; //Assume we have an array with all zeros                                                                                                                                                        
         bool first20hasZeros = false; //Checking header                                                                                                                                                                       

         if(N ==0){ //If payloadwords = 0 then the array is empty                                                                
           ++emptyZS_arrays;
           std::cout<< "Found empty ZS array at ii = "<<ii<< " \n";// extra check ->Raw has no empty arrays it seems                                                                                                          
         }else{
           for(size_t k = 0; k< N; ++k){
             if(payloadarray[k] != 0){
               allZeros = false; //Break if you find an element that isn't zero                                                                       
               break;
             }
           }
           //Check to see if first 20 elements have any zero                                                                                                                                                                   
           for(size_t z20 =0; z20 < std::min<size_t>(N,20); ++z20){
             if(payloadarray[z20] == 0){
               first20hasZeros = true;//Break if you find any element with a zero, will disrupt header                                                                                                                         
               break;
             }
           }

           if (allZeros == true || first20hasZeros == true){
             ++zeroZS_arrays; //Adds to the zero array counter                                                                                                                                                               
             std::cout<<"PayloadWords for zeroZS array: "<<N<< std::endl; //extra check                                                                                                                                       
             for(size_t zz = 0; zz < std::min<size_t>(N,10);++zz){
               std::cout<<payloadarray[zz]<<" "; //Prints out first 10 elemenmts to check we have all zeros                                                                                                                    
             }
             std::cout<<" ZS ii = " << ii<< "\n";

           } else {
             for(size_t kk = 0; kk < std::min<size_t>(N,20); ++kk){
               std::cout<<payloadarray[kk]<<" "; // prints out the first 20 elements                                                                                                                                           
             }
             std::cout<<" ZS ii = "<< ii <<"\n";
             stm_waveform.set_data(stm_frag.payloadWords(), stm_frag.payloadBegin());
             zs_waveform_digis->emplace_back(stm_waveform);
           }
         }

	 //end of check         
  
	 //stm_waveform.set_data(stm_frag.payloadWords(), stm_frag.payloadBegin());
	 //zs_waveform_digis->emplace_back(stm_waveform);
       
       }
       else if (stm_frag.isMWD()) {
	 int n_MWD_digis = stm_frag.payloadWords(); //number read -> to digis
	 for (int i_MWD = 0 ; i_MWD < n_MWD_digis; ++i_MWD){
	   auto const* pointer  = stm_frag.payloadBegin(); //tells where to read data
	   int16_t i_pointer = pointer[i_MWD]; //Retrives value of ith index of pointer
	   mu2e::STMMWDDigi mwd_digi(0,i_pointer);
	   mwd_digis->emplace_back(mwd_digi);
	   std::cout<< "mwd_digis size = " << mwd_digis->size()<<"\n";
	 }

       }
    }
  }

  //Final checks
  std::cout<<"---------------------------------------------------------\n"<<"\n";
  std::cout << "raw_waveform_digis total size = "<<raw_waveform_digis->size()<<"\n"; //Check before being put in
  std::cout << "zs_waveform_digis total size = "<<zs_waveform_digis->size()<<"\n";
  std::cout << "mwd_digis total size = " << mwd_digis->size()<<"\n";

  std::cout<<"---------------------------------------------------------\n"<<"\n";
  std::cout<< "Number of Zero arrays for Raw = " <<zeroRaw_arrays<<"\n";
  std::cout<< "Number of Zero arrays for ZS = " <<zeroZS_arrays<<"\n";

  std::cout<<"---------------------------------------------------------\n"<<"\n";
  std::cout<< "Number of empty arrays for Raw = " <<emptyRaw_arrays<<"\n";
  std::cout<< "Number of empty arrays for ZS = " <<emptyZS_arrays<<"\n";


  
  event.put(std::move(raw_waveform_digis), "raw");
  event.put(std::move(zs_waveform_digis), "zs");
  event.put(std::move(mwd_digis), "mwd");
    
} // produce()

// ======================================================================

DEFINE_ART_MODULE(STMDigisFromFragments)
