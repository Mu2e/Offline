#ifndef MCDataProducts_CaloShowerStep_hh
#define MCDataProducts_CaloShowerStep_hh

#include "MCDataProducts/inc/SimParticle.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include <iostream>


namespace mu2e {


   class CaloShowerStep 
   {
       public:	 
	  CaloShowerStep(): sim_(),volumeId_(-1),nCompress_(0),time_(0),energyDepG4_(0),energyDepBirks_(0),momentumIn_(0),pos_() {}
	  
	  CaloShowerStep(int volumeId,      const art::Ptr<SimParticle>& sim, int nCompress,    float time,  
	                 float energyDepG4, float energyDepBirks,             float momentumIn, const CLHEP::Hep3Vector& pos): 
	    sim_(sim),
	    volumeId_(volumeId),
	    nCompress_(nCompress),
	    time_(time),
	    energyDepG4_(energyDepG4),
	    energyDepBirks_(energyDepBirks),
	    momentumIn_(momentumIn),
	    pos_(pos)
	  {}

	  const art::Ptr<SimParticle>&   simParticle()    const {return sim_;}
	  const CLHEP::Hep3Vector&       position()       const {return pos_;}	  
	        int                      volumeId()       const {return volumeId_;}
	        int                      nCompress()      const {return nCompress_;}
	        float                    time()           const {return time_;}
	        float                    energyDepG4()    const {return energyDepG4_;}
	        float                    energyDepBirks() const {return energyDepBirks_;}
	        float                    momentumIn()     const {return momentumIn_;}

          void setSimParticle(const art::Ptr<SimParticle>& sim) {sim_ = sim;}                                                    
          
          void print(std::ostream& ost = std::cout) const {ost<<"Calo Shower content    volumeId = "<<volumeId_<<"  pid="<<sim_->pdgId()
	                                                      <<"  edep="<<energyDepG4_<<"  evis="<<energyDepBirks_
                                                              <<"  time="<<time_       <<"  position = "<<pos_<<std::endl;}
       private:       
	    art::Ptr<SimParticle>   sim_;            
	    int                     volumeId_;       
	    int                     nCompress_;      
	    float                   time_;           
	    float                   energyDepG4_;    
	    float                   energyDepBirks_; 
	    float                   momentumIn_;     
	    CLHEP::Hep3Vector       pos_;            
   };



   typedef std::vector<mu2e::CaloShowerStep> CaloShowerStepCollection;

} 

#endif
