//
//
// Author: Lisa Goodeough
// Date: 2017/07/09
//
//

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4WorkerInitialization.hh"

//art includes
#include "fhiclcpp/ParameterSet.h"

using namespace std;

namespace mu2e {

Mu2eG4WorkerInitialization::Mu2eG4WorkerInitialization(const fhicl::ParameterSet& pset
                                 )
    :
    G4UserWorkerInitialization(),
    pset_(pset)
    {
        
        std::cout << "We are in the c'tor of WorkerInitialize!!!" << std::endl;
        
    }

Mu2eG4WorkerInitialization::~Mu2eG4WorkerInitialization()
    {}
    
    
void Mu2eG4WorkerInitialization::WorkerInitialize() const
    {
        std::cout << "We are at WorkerInitialize!!!" << std::endl;
    }
    
void Mu2eG4WorkerInitialization::WorkerStart() const
    {
        std::cout << "We are at WorkerStart!!!" << std::endl;
    }

void Mu2eG4WorkerInitialization::WorkerRunStart() const
    {
        std::cout << "We are at WorkerRunStart!!!" << std::endl;
    }

void Mu2eG4WorkerInitialization::WorkerRunEnd() const
    {
        std::cout << "We are at WorkerRunEnd!!!" << std::endl;
    }
    
void Mu2eG4WorkerInitialization::WorkerStop() const
    {
        std::cout << "We are at WorkerStop!!!" << std::endl;
    }
}  // end namespace mu2e



