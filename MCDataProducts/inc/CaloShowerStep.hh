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
          CaloShowerStep(): sim_(),volumeG4ID_(-1),nCompress_(0),time_(0),energyDepG4_(0),energyDepBirks_(0),momentumIn_(0),pos_() {}
          
          CaloShowerStep(int volumeG4ID, const art::Ptr<SimParticle>& sim, unsigned nCompress, float time, float energyDepG4, 
                         float energyDepBirks, float momentumIn, const CLHEP::Hep3Vector& pos): 
            sim_(sim),
            volumeG4ID_(volumeG4ID),
            nCompress_(nCompress),
            time_(time),
            energyDepG4_(energyDepG4),
            energyDepBirks_(energyDepBirks),
            momentumIn_(momentumIn),
            pos_(pos)
          {}

          const art::Ptr<SimParticle>&   simParticle()    const {return sim_;}
          const CLHEP::Hep3Vector&       position()       const {return pos_;}          
          int                            volumeG4ID()     const {return volumeG4ID_;}
          unsigned                       nCompress()      const {return nCompress_;}
          double                         time()           const {return time_;}
          float                          energyDepG4()    const {return energyDepG4_;}
          float                          energyDepBirks() const {return energyDepBirks_;}
          float                          momentumIn()     const {return momentumIn_;}

          void setSimParticle(const art::Ptr<SimParticle>& sim) {sim_ = sim;}                                                    
          
          void print(std::ostream& ost = std::cout) const {ost<<"Calo Shower content    volumeG4ID = "<<volumeG4ID_<<"  pid="<<sim_->pdgId()
                                                              <<"  edep="<<energyDepG4_<<"  evis="<<energyDepBirks_
                                                              <<"  time="<<time_       <<"  position = "<<pos_<<std::endl;}
       private:       
          art::Ptr<SimParticle>   sim_;            
          int                     volumeG4ID_;       
          unsigned                nCompress_;      
          double                  time_;           
          float                   energyDepG4_;    
          float                   energyDepBirks_; 
          float                   momentumIn_;     
          CLHEP::Hep3Vector       pos_;            
   };


   using CaloShowerStepCollection = std::vector<mu2e::CaloShowerStep>;
} 

#endif
