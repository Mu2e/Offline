//
//
// Author: Lisa Goodeough
// Date: 2017/08/07
//
//

#include "Mu2eG4/inc/Mu2eG4MasterRunAction.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

using namespace std;

namespace mu2e {

Mu2eG4MasterRunAction::Mu2eG4MasterRunAction(PhysicalVolumeHelper* phys_volume_helper,
                                             //PhysicsProcessInfo* phys_process_info
                                             std::vector< PhysicsProcessInfo >* phys_process_info_vec)
    :
    G4UserRunAction(),
    _physVolHelper(phys_volume_helper),
    //_processInfo(phys_process_info)
    PhysicsProcessInfoVector(phys_process_info_vec)
    {
        std::cout << "AT Mu2eG4MasterRunAction c'tor" << std::endl;
    }
    
    
Mu2eG4MasterRunAction::~Mu2eG4MasterRunAction()
    {
        std::cout << "AT Mu2eG4MasterRunAction destructor" << std::endl;
    }
    
    
void Mu2eG4MasterRunAction::BeginOfRunAction(const G4Run* aRun)
    {
        //this class is ONLY called in MT mode
        //we want these actions performed only in the Master thread
        _physVolHelper->beginRun();//map w/~20,000 entries
        
        for (unsigned i = 0; i < PhysicsProcessInfoVector->size(); i++) {
            PhysicsProcessInfoVector->at(i).beginRun();
            
            //std::cout << "calling PPI.beginRun() for PPI.at " << i << std::endl;
            
        }
        //_processInfo->beginRun();

    }
    
    
void Mu2eG4MasterRunAction::EndOfRunAction(const G4Run* aRun)
    {
        std::cout << "calling Mu2eG4MasterRA::EndOfRunAction()" << std::endl;
        
        for (unsigned i = 0; i < PhysicsProcessInfoVector->size(); i++) {
            PhysicsProcessInfoVector->at(i).endRun();
        }
        //_processInfo->endRun();
        
        //nothing for this
        // _physVolHelper.endRun();
        
    }

    

}  // end namespace mu2e



