// The method SimParticle::parent() allows to navigate the provenance
// chain up.  However by just following the parent() method one only
// gets particles produced in one jobs.  If a simulation job was
// factorized so that SimParticles were tracked to a VirtualDetector
// and stopped, then another G4 job used the FromStepPointMCs module
// to read in VD hits and continue simulation, the full provenance chain
// looks like
//
//  GenParticle -> SimParticle [...] -> StepPointMC  <=(Assn)=>  GenParticle2 -> SimParticle [...]
//
// where (Assn) is the art::Association object used to connect
// GenParticle2 in the second job to the previous history.
//
// The SimParticleParentGetter class allows transparent
// navigation through the chain of SimParticles up to the
// most primary SimParticle regardless of whether and how the
// simulation job was split in pieces and restarted.
//
// Usage:
//
//    // once per event, outside of loops
//    SimParticleParentGetter pg(event);
//    .....
//    // In a loop or anywhere
//    art::Ptr<SimParticle> current = ...;
//    art::Ptr<SimParticle> parent = pg.parent(current);
//
// Andrei Gaponenko, 2012
// The class has been extended to get the the parent of a GenParticle
// that comes form the end point of a previous job SimParticle
// created with the modue EventGenerator/src/FromSimParticleEndPoint_module.cc



#ifndef Mu2eUtilities_SimParticleParentGetter_hh
#define Mu2eUtilities_SimParticleParentGetter_hh

#include <map>

#include "canvas/Persistency/Common/Ptr.h"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StepPointMC.hh"

namespace art { class Event; }

namespace mu2e {

  class SimParticleParentGetter {
    const art::Event *evt_;

    typedef std::map<art::Ptr<GenParticle>, art::Ptr<StepPointMC> > StepPointMapType;
    mutable StepPointMapType stepPointMap_;

    typedef std::map<art::Ptr<GenParticle>, art::Ptr<SimParticle> > SimParticleMapType;
    mutable SimParticleMapType simParticleMap_;

    static art::Ptr<SimParticle> nullPtr_;
  public:
    SimParticleParentGetter(const art::Event& evt);

    // Returns null when the end of the chain is reached.
    const art::Ptr<SimParticle>& parent(const art::Ptr<SimParticle>& particle) const;
  };
}

#endif/*Mu2eUtilities_SimParticleParentGetter_hh*/
