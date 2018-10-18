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
                                             std::vector< PhysicsProcessInfo >* phys_process_info_vec)
    :
    G4UserRunAction(),
    _physVolHelper(phys_volume_helper),
    PhysicsProcessInfoVector(phys_process_info_vec)
    {}
    
    
Mu2eG4MasterRunAction::~Mu2eG4MasterRunAction()
    {}
    
    
void Mu2eG4MasterRunAction::BeginOfRunAction(const G4Run* aRun)
    {
        //this class is ONLY called in MT mode
        //we want these actions performed only in the Master thread
        _physVolHelper->beginRun();//map w/~20,000 entries
        
        for (unsigned i = 0; i < PhysicsProcessInfoVector->size(); i++) {
            PhysicsProcessInfoVector->at(i).beginRun();
        }
    }
    
    
void Mu2eG4MasterRunAction::EndOfRunAction(const G4Run* aRun)
    {        
        for (unsigned i = 0; i < PhysicsProcessInfoVector->size(); i++) {
            PhysicsProcessInfoVector->at(i).endRun();
        }
    }

    

}  // end namespace mu2e



