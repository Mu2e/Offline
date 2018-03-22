//
// Struct describing a single straw hit assigned to a track, for use in TTree diagnostics
// $Id: TrkStrawHitInfo.hh,v 1.2 2014/09/22 12:13:17 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/09/22 12:13:17 $
//
#ifndef TrkStrawHitInfo_HH
#define TrkStrawHitInfo_HH
#include "Rtypes.h"
#include "CLHEP/Vector/ThreeVector.h"
namespace mu2e 
{
  struct TrkStrawHitInfo {
    Bool_t _active;   // is this hit used in the track?
    Bool_t _dhit;     // is this hit one of a pair in a panel?
    Bool_t _dactive;  // is this hit one of a panel pair that are both used in the track
    Int_t _plane, _panel, _layer, _straw;  // StrawId fields for the straw hit
    Int_t _ambig;     // left-right amibiguity.  This signes the angular momentum of the track WRT the wire
    Int_t _driftend; // which end(s) was used in computing the drift
    Float_t _tdcal, _tdhv; // drift times for the 2 ends
    Float_t _totcal, _tothv; // TOT for the 2 ends
    CLHEP::Hep3Vector _poca; // Position of Point Of Closest Approach (POCA)
    Float_t _resid;	  // residual = Distance Of Closest Approach (DOCA) between the drift cylinder and the track, signed by the track angular momentum WRT the wire
    Float_t _residerr;	  // error on the residual, including components from the hit, the track, and potentially other effects 
    Float_t _rdrift, _rdrifterr;	  // drift radius and error of this hit
    Float_t _trklen;	// distance along the helix of the POCA for this hit
    Float_t _doca;	// DOCA of this hit
    Float_t _exerr, _penerr; // external (simulated annealing) and penalty (pat. rec.) errors added to this hit in the total residual error
    Float_t _t0, _t0err;  // Estimated time the particle passed the POCA of this hit (and error).  Note this is mass hypothesis dependent
    Float_t _ht;    // reconstructed time (including drift time) of this hit
    Float_t _dt, _tddist, _tdderr;	// Delta-t, Distance from the straw center predicted by the time difference measurement (and error)
    Float_t _hlen;    // length along the straw ffrom the straw center of the POCA
    Float_t _wdot;    // dot-product of the track direction at POCA with the wire direction
    Float_t _bdot;    // dot product of the track direction at POCA with the local BField
    Float_t _edep;    // reconstructed energy deposition from ADC measurement
    Float_t _dx;      // estimated distance through the straw gas traversed by this particle, given the DOCA and track parameters
    TrkStrawHitInfo() : _active(false), _dhit(false), _dactive(false), 
    _plane(-1), _panel(-1), _layer(-1), _straw(-1), _ambig(-1),_driftend(-1),
    _tdcal(-1000.0),_tdhv(-1000.0),
    _resid(-1000.0), _residerr(-1000.0), _rdrift(-1000.0), _rdrifterr(-1000.0),_trklen(-1000.0),
    _doca(-1000.0), _exerr(-1000.0), _penerr(-1000.0), _t0(-1000.0), _t0err(-1000.0),
    _ht(-1000.0), _dt(-1000.0), _tddist(-1000.0), _tdderr(-1000.0), _hlen(-1000.0), _wdot(-1000.0), _bdot(-1000.0),
    _edep(-1000.0), _dx(-1000.0)  {}
  };
}
#endif
