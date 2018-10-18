//
// An EDAnalyzer module that verifies that the GenParticle to SimParticle bookkeeping has been done correctly in the G4MT version of the code
//
// Original author Lisa Goodenough
//



#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
//#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include <iostream>
#include <string>
#include <fstream>


using namespace std;


namespace mu2e {

  class MTVerification : public art::EDAnalyzer {
  public:

      explicit MTVerification(fhicl::ParameterSet const& pset);
      virtual ~MTVerification() { }

      virtual void beginJob() override;
      virtual void endJob() override;
      
      //virtual void beginRun(art::Run &r);
      
      virtual void beginRun(art::Run const& r ) override;

      virtual void analyze(const art::Event& e) override;

  private:
 
 
      // Input tag of the generated particles.
      art::InputTag _gensTag;

      // Module label of the g4 module that made the hits.
      string _g4ModuleLabel;

      // Module label of the generator module that was passed as input to G4.
      string _generatorModuleLabel;
      
      ofstream myOutFile;
      
    // Examine GenParticle and SimParticle Collections for consistency
    void checkGenParticleInSimParticle ( const art::Event& event );
      

  };
    
  MTVerification::MTVerification(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    // Run time parameters
    _gensTag(pset.get<string>("gensTag","generate")),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _generatorModuleLabel(pset.get<string>("generatorModuleLabel"))
    
    {
      }

    
  void MTVerification::beginJob(){
      
      myOutFile.open("/home/goodenou/Mu2e_MT/Offline/MTVerification.txt");
  }
    
    
  void MTVerification::beginRun(art::Run const& run){
        
        myOutFile << "runID: " << run.id() << endl;
    }


  void MTVerification::analyze(const art::Event& event) {
 
      // Call code that analyzes the Sim and Gen Particle Collections
       checkGenParticleInSimParticle(event);

  }
  
  void MTVerification::checkGenParticleInSimParticle(const art::Event& event){
      
      bool areEqual = true;
      
      myOutFile  << "******************************************" << endl;
      
      // Get a ValidHandle to the generated particles.
      //auto gens = event.getValidHandle<GenParticleCollection>(_gensTag);
      art::ValidHandle<GenParticleCollection> genParticles = event.getValidHandle<GenParticleCollection>(_gensTag);
 
      art::ValidHandle<SimParticleCollection> simParticles = event.getValidHandle<SimParticleCollection>(_g4ModuleLabel);
      
      //art::Handle<SimParticleCollection> simParticles;
      //event.getByLabel(_g4ModuleLabel,simParticles);
      //bool haveSimPart = ( simParticles.isValid() );

      for ( auto const& gen : *genParticles ){
          
          myOutFile  << "MTVerification Gen: evtID="
                << event.id().event() << ", "
                << gen.pdgId() << " "
                << gen.generatorId() << " "
                << gen.time() << " "
                << gen.properTime() << " "
                << gen.position() << " "
                << gen.momentum() << " "
                << endl;
      }


      for ( SimParticleCollection::const_iterator i=simParticles->begin();
            i!=simParticles->end(); ++i ){
          
          if (i->first.asInt() == 1 ) {
              
              SimParticle const& sim = i->second;
              
              // Information about generated particle.
              GenParticle const& gen = *sim.genParticle();
              GenId genId(gen.generatorId());
              
              myOutFile  << "MTVerification Sim: evtID="
                    << event.id().event()       << ", "
                    << sim.pdgId()              << " "
                    << genId                    << " "
                    << sim.startGlobalTime()    << " "
                    << sim.startProperTime()    << " "
                    << sim.startPosition()      << " "
                    << sim.startMomentum()
                    << endl;
              
              if ( sim.pdgId() != gen.pdgId() ){
                  areEqual = false;
                  myOutFile << "pdgId not the same" << endl;
              }
              
              if ( genId != gen.generatorId() ){
                  areEqual = false;
                  myOutFile << "genId not the same" << endl;
              }
              
              if ( sim.startGlobalTime() != gen.time() ){
                  areEqual = false;
                  myOutFile << "global time not the same" << endl;
              }
              
              if ( sim.startProperTime() != gen.properTime() ){
                  areEqual = false;
                  myOutFile << "proper time not the same" << endl;
              }
              
              if ( sim.startPosition() != gen.position() ){
                  areEqual = false;
                  myOutFile << "position not the same" << endl;
              }
              
              if ( sim.startMomentum() != gen.momentum() ){
                  areEqual = false;
                  myOutFile << "momentum not the same" << endl;
              }
              
              if (areEqual) {
                  myOutFile  << "GenParticle and 1st SimParticle are the same!" << endl;
              } else {
                  myOutFile  << "GenParticle and 1st SimParticle are NOT the same!" << endl;
              }
              
          }

      }
      
      
      

  } // end checkGenParticleInSimParticle


  void MTVerification::endJob(){
      
      //myOutFile.close();

  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::MTVerification);

