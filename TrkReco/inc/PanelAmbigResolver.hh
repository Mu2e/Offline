//
// class to resolve hit ambiguities by panel, assuming a reasonable track
// fit as input
//
// Original author: David Brown (LBNL), 2012
//
//
#ifndef mu2e_PanelAmbig_PanelAmbigResolver_HH
#define mu2e_PanelAmbig_PanelAmbigResolver_HH
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/AmbigResolver.hh"
#include "TrkReco/inc/PanelAmbigStructs.hh"
#include "DataProducts/inc/PanelId.hh"
#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/
#include <cstddef>
#include <vector>
// Root
#include "Rtypes.h"
#include "TTree.h"

class KalRep;
class TrkT0;

namespace mu2e {
  namespace PanelAmbig {
    class PanelAmbigResolver : public AmbigResolver {
      public:
	// construct from parameter set
#ifndef __GCCXML__
	explicit PanelAmbigResolver(fhicl::ParameterSet const&, double tmpErr, size_t iter );
#endif/*__GCCXML__*/
	virtual ~PanelAmbigResolver();
	// resolve a track.  Depending on the configuration, this might
	// update the hit state and the t0 value.
	virtual bool resolveTrk(KalRep* krep) const;
      private:
	// resolve the ambiguity on a single panel
	bool resolvePanel(TrkStrawHitVector& phits, KalRep* krep) const;
	// fill information about a given panel's track and hits
	bool fillPanelInfo(TrkStrawHitVector const& phits, const KalRep* krep, PanelInfo& pinfo) const;
	// compute the panel result for a given ambiguity/activity state and the ionput t0
	void fillResult(PanelInfo const& pinfo,TrkT0 const& t0, PanelResult& result) const;
	// parameters
	double _minsep; // minimum chisquared separation between best solution and the rest to consider a panel resolved
	double _inactivepenalty; // chisquared penalty for an inactive hit
	double _nullpenalty; // chisquared penalty for a null ambiguity hit
	double _penaltyres; // resolution term to add to hits if ambiguity/activity can't be resolved
	double _trkpenaltyres; // resolution term to add to track
	double _nullerr2; // additional error (squared) for hits with null ambiguity
	bool _addtrkpos; // add constraint associated with the track position to the chisquared
	HSV _allowed; // allowed states of a TrkStrawHit
	double _maxhitu; // maximum u value allowed for a hit
	bool _fixunallowed; // fix the state of any hit whose initial state isn't allowed
	unsigned _maxnpanel; // max # of hits to consider for a panel
	int _diag; // diagnostic level`
	// TTree variables, mutable so they don't change const
	mutable TTree *_padiag, *_pudiag; // diagnostic TTree
	mutable Int_t _nuhits, _nrhits; // # hits in this panel
	mutable Int_t _nactive; // # active hits in this panel
	mutable Int_t _nused; // # hits used to compute chisquared
	mutable Int_t _nres; // # of results for this panel (combinatorics)
	mutable Float_t _tupos; // u track position
	mutable Float_t _tuerr; // u track position
	mutable TSHUIV _uinfo; // u position of hits in panel
	mutable Float_t _mctupos;
	mutable PRV _results;
    };
  } // PanelAmbig namespace
} // mu2e namespace

#endif

