#ifndef MCDataProducts_CrvStep_hh
#define MCDataProducts_CrvStep_hh
//
// Class to summarize the passage of a single particle through a single CRV counter
// This consolidates the G4 steps and insulates the downstream CRV response simulation 
// from details of the G4 model
//
#include "canvas/Persistency/Common/Ptr.h"

#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <Rtypes.h>

namespace mu2e 
{
  class CrvStep 
  {
    public:

      CrvStep() : _visibleEDep(0.0), _startTime(0.0), _endTime(0.0), _pathLength(0.0) {}
      CrvStep(CRSScintillatorBarIndex barIndex, float visibleEDep, double startTime, double endTime,
              const XYZVec &startPos, const XYZVec &endPos, const XYZVec &startMom, float endMom, 
              float pathLength, art::Ptr<SimParticle> const& simParticle) :
              _barIndex(barIndex), _visibleEDep(visibleEDep), 
              _startTime(startTime), _endTime(endTime), _startPos(startPos), _endPos(endPos),
              _startMom(startMom), _endMom(endMom),
              _pathLength(pathLength), _simParticle(simParticle) {}

      CRSScintillatorBarIndex      barIndex() const    {return _barIndex;}
      float                        visibleEDep() const {return _visibleEDep;}
      double                       startTime() const   {return _startTime;}
      double                       endTime() const     {return _endTime;}

      XYZVec const&                startPos() const    {return _startPos;}
      XYZVec const&                endPos() const      {return _endPos;}
      XYZVec const&                startMom() const    {return _startMom;}
      float                        endMom() const      {return _endMom;}

      CLHEP::Hep3Vector            startPosition() const {return Geom::Hep3Vec(_startPos);}
      CLHEP::Hep3Vector            endPosition() const   {return Geom::Hep3Vec(_endPos);}
      CLHEP::Hep3Vector            startMomentum() const {return Geom::Hep3Vec(_startMom);}

      float                        pathLength() const  {return _pathLength;}

      art::Ptr<SimParticle> const& simParticle() const {return _simParticle;}
      art::Ptr<SimParticle>&       simParticle()       {return _simParticle;}

    private:
      CRSScintillatorBarIndex _barIndex;
      float                   _visibleEDep; 
      double                  _startTime, _endTime; //must be double to allow for long-lived particles
      XYZVec                  _startPos,  _endPos;
      XYZVec                  _startMom;
      float                   _endMom;              //direction of end momentum is not known
      float                   _pathLength;          //the actual step length, which may be longer 
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

