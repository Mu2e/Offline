//
// BaBar hit object corresponding to a single cell hit
//
// $Id: TrkCellHit.hh,v 1.2 2012/12/04 00:51:26 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:26 $
//
// Original author G. Tassielli
//
#ifndef TrkCellHit_HH
#define TrkCellHit_HH

#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTestsI/inc/DetCellGasElem.hh"
#include "KalmanTestsI/inc/DetFieldWireElem.hh"
#include "KalmanTestsI/inc/DetSenseWireElem.hh"
//#include "KalmanTests/inc/DetStrawHitType.hh"

#include "ITrackerGeom/inc/ITracker.hh"

#include <string>

namespace mu2e 
{
  class TrkCellHit : public TrkStrawHit {

          friend class DetFieldWireElem;
          friend class DetSenseWireElem;

  public:

    TrkCellHit(const StrawHit& strawhit, const Straw& straw,unsigned istraw,
    const TrkT0& trkt0, double fltlen, double exterr, double maxdriftpull);
    TrkCellHit(const StrawHit& strawhit, const Straw& straw,unsigned istraw,
    const TrkT0& trkt0, double fltlen, double exterr, double maxdriftpull, std::string matOnly);
    //TrkCellHit(const TrkCellHit& other, const TrkDifTraj* trkTraj);
    virtual ~TrkCellHit();
//  Simplistic implementation of TrkHitOnTrk interface.  Lie where necessary
    virtual TrkCellHit* clone(TrkRep* parentRep, const TrkDifTraj* trkTraj = 0) const;
    const double cellPath(double pftl,Hep3Vector const& tdir, HepPoint const& tpos, double omega=0.0, double cosDip=1.0, double angle=0.0); // track pathlength through the field wires of the cell
    const double senseWirePath(double pftl,Hep3Vector const& tdir, HepPoint const& tpos, double omega=0.0, double cosDip=1.0, double angle=0.0); // track pathlength through the field wires of the cell
    const double fieldWirePath(double pftl,Hep3Vector const& tdir, HepPoint const& tpos, double omega=0.0, double cosDip=1.0, double angle=0.0); // track pathlength through the field wires of the cell
    const double cellGasPath(double pftl,Hep3Vector const& tdir, HepPoint const& tpos, double omega=0.0, double cosDip=1.0, double angle=0.0); // track pathlength through the gas of the cell
    double entryDeltaPath()    const { return _entryDltPath;   }
    double exitDeltaPath()     const { return _exitDltPath;    }
    double centerDeltaPath()   const { return _centerDlPath;   }
    bool   isGrowingRad()      const { return _growingRad;     }
    bool   pcaOutOfCut()       const { return _hitOutOfCut;    }
    bool   hitSwire()          const { return _hitSWire;       }
    bool   hitFwireTop()       const { return _hitFWire_t;     }
    bool   hitFwireSide()      const { return _hitFWire_s;     }
    bool   hitFwireBot()       const { return _hitFWire_b;     }
    double getSwirePath()      const { return _pathlInSWire;   }
    double getFwirePathTop()   const { return _pathlInFWire_t; }
    double getFwirePathSide()  const { return _pathlInFWire_s; }
    double getFwirePathBot()   const { return _pathlInFWire_b; }
    int    GetSuperLayer() const;
    int    GetCelRing() const;
    int    GetCell() const;

// access to associated detector elements
    DetCellGasElem const& cellGasElem() const { return _cgelem; }
    DetFieldWireElem const& fwireElemBottom() const { return _fwelem_b; }
    DetFieldWireElem const& fwireElemSide() const { return _fwelem_s; }
    DetFieldWireElem const& fwireElemTop() const { return _fwelem_t; }
    DetSenseWireElem const& swireElem() const { return _swelem; }
    CellGeometryHandle *cellHandle() const;
  protected:
    TrkCellHit(const TrkCellHit& other, TrkRep* rep);
    virtual void updateSignalTime();
    virtual void updateDrift();

// DetModel stuff
    static DetStrawHitType* cgtype(std::string gasMat="ITgasAuto");
    static DetStrawHitType* fwtype(std::string fwMat="ITFWireAuto");
    static DetStrawHitType* swtype(std::string swMat="ITSWireAuto");
    DetCellGasElem _cgelem;
    DetFieldWireElem _fwelem_b;
    DetFieldWireElem _fwelem_s;
    DetFieldWireElem _fwelem_t;
    DetSenseWireElem _swelem;
    CellGeometryHandle *_itwp;

    bool _evalCellPath;
    double _entryDltPath;
    double _exitDltPath;
    double _centerDlPath;
    bool _growingRad;
    bool _hitOutOfCut;

    bool   _hitSWire;
    bool   _hitFWire_t;
    bool   _hitFWire_s;
    bool   _hitFWire_b;
    double _pathlInSWire;
    double _pathlInFWire_t;
    double _pathlInFWire_s;
    double _pathlInFWire_b;
    //HelixTraj *_hTdef;

  };
//// unary functor to select TrkCellHit from a given hit
//  struct FindTrkCellHit {
//    FindTrkCellHit(StrawHit const& strawhit) : _strawhit(strawhit) {}
//    bool operator () (TrkCellHit* const& tsh ) { return tsh->strawHit() == _strawhit; }
//    StrawHit const& _strawhit;
//  };

  
}

#endif
