#ifndef MCDataProducts_CrvStep_hh
#define MCDataProducts_CrvStep_hh
//
// Class to summarize the passage of a single particle through a single CRV counter
// This consolidates the G4 steps and insulates the downstream CRV response simulation
// from details of the G4 model
//
#include "canvas/Persistency/Common/Ptr.h"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <Rtypes.h>

namespace mu2e
{
  class CrvStep
  {
    public:

      CrvStep() {}
      CrvStep(CRSScintillatorBarIndex barIndex, float visibleEDep, double startTime, double endTime,
              const XYZVectorF &startPos, const XYZVectorF &endPos, const XYZVectorF &startMom, float endMom,
              float pathLength, art::Ptr<SimParticle> const& simParticle) :
              _barIndex(barIndex), _visibleEDep(visibleEDep),
              _startTime(startTime), _endTime(endTime), _startPos(startPos), _endPos(endPos),
              _startMom(startMom), _endMom(endMom),
              _pathLength(pathLength), _simParticle(simParticle) {}

      CRSScintillatorBarIndex      barIndex() const    {return _barIndex;}
      float                        visibleEDep() const {return _visibleEDep;}
      double                       startTime() const   {return _startTime;}
      double                       endTime() const     {return _endTime;}
      double&                       startTime() {return _startTime;}; // non-const used for resampling
      double&                       endTime() {return _endTime;}

      XYZVectorF const&                startPos() const    {return _startPos;}
      XYZVectorF const&                endPos() const      {return _endPos;}
      XYZVectorF const&                startMom() const    {return _startMom;}
      float                        endMom() const      {return _endMom;}

      CLHEP::Hep3Vector            startPosition() const {return GenVector::Hep3Vec(_startPos);}
      CLHEP::Hep3Vector            endPosition() const   {return GenVector::Hep3Vec(_endPos);}
      CLHEP::Hep3Vector            startMomentum() const {return GenVector::Hep3Vec(_startMom);}

      float                        pathLength() const  {return _pathLength;}

      art::Ptr<SimParticle> const& simParticle() const {return _simParticle;}
      art::Ptr<SimParticle>&       simParticle()       {return _simParticle;}

    private:
      CRSScintillatorBarIndex _barIndex;
      float                   _visibleEDep{0};
      double                  _startTime{0}, _endTime{0}; //must be double to allow for long-lived particles
      XYZVectorF              _startPos,  _endPos;
      XYZVectorF              _startMom;
      float                   _endMom{0};              //direction of end momentum is not known
      float                   _pathLength{0};          //the actual step length, which may be longer
                                                    //than the differences between endPos and startPos
      art::Ptr<SimParticle>   _simParticle;
  };

  typedef std::vector<CrvStep> CrvStepCollection;

  inline std::ostream& operator<<( std::ostream& ost, CrvStep const& cs){
    ost << "CrvStep BarIndex " << cs.barIndex()
    << " visible energy deposit " << cs.visibleEDep()
    << " path length " << cs.pathLength();
    return ost;
  }

}

#endif

