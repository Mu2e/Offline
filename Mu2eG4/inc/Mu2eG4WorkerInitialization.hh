#ifndef Mu2eG4_WorkerInit_hh
#define Mu2eG4_WorkerInit_hh
//
// Mu2eG4WorkerInitialization.hh provides declarations for the built-in UserWorkerInitialization
// class for the Mu2e G4 simulation.
//
// Author: Lisa Goodeough
// Date: 2019/07/09
//
//


//G4 includes
#include "G4UserWorkerInitialization.hh"

//C++ includes
#include <iostream>

namespace fhicl { class ParameterSet; }

namespace mu2e {

class Mu2eG4WorkerInitialization : public G4UserWorkerInitialization
{
  public:
    Mu2eG4WorkerInitialization(const fhicl::ParameterSet& pset
                    );

    virtual ~Mu2eG4WorkerInitialization();

    //methods
    virtual void WorkerInitialize() const;
    virtual void WorkerStart() const;
    virtual void WorkerRunStart() const;
    virtual void WorkerRunEnd() const;
    virtual void WorkerStop() const;

  private:
    //data members

    const fhicl::ParameterSet& pset_;
    //const bool use_G4MT_;

};

}  // end namespace mu2e
#endif /* Mu2eG4_WorkerInit_hh */
