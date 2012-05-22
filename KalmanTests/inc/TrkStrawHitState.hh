//
// Object to allow exhaustively iterating over state permutations of a set of TrkStrawHits
//
// $Id: TrkStrawHitState.hh,v 1.2 2012/05/22 21:35:42 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/22 21:35:42 $
//
#ifndef TrkStrawHitState_HH
#define TrkStrawHitState_HH
#include <cstddef>
#include <vector>

namespace mu2e {
  class TrkStrawHit;
// class to describe the ambiguity and activity state of a single TrkStrawHit
  class TrkStrawHitState {
    public:
      enum TSHState {ignore=-1,noambig=0,negambig,posambig,inactive,nstates};
      // set state explicitly.  
      TrkStrawHitState(TSHState state=ignore, TrkStrawHit* tsh=0) : _state(state), _tsh(tsh) {}
      // set state from an existing strawhit
      TrkStrawHitState(TrkStrawHit* tsh);
      // acces the underlying hit
      TrkStrawHit* hit() { return _tsh; }
      const TrkStrawHit* hit() const { return _tsh; }
      // set the state of the straw hit to correspond to its cached hit
      void setHitState() const;
      // return the state as an integer
      TSHState state() const { return _state; }
      // total number of states
      static unsigned nStates() { return nstates; }
      // return to default state
      void setState(TSHState state) { _state = state; }
      // comparators
      bool operator == (TrkStrawHitState const& other) const { return _state == other._state; }
      bool operator != (TrkStrawHitState const& other) const { return !operator==(other); }
    private:
      TSHState _state;
      TrkStrawHit* _tsh; // cache of strawhit, used to validate
  };

  typedef std::vector<TrkStrawHit*> TSHV;
  typedef std::vector<TrkStrawHitState> TSHSV;
  typedef std::vector<TrkStrawHitState::TSHState> TSHSSV;

  // class to describe the states of a vector of TrkSTrawHits
  class TrkStrawHitStateVector {
    public:
      TrkStrawHitStateVector();  // default constructor allows NO STATES!!!
      // construct from a vector of TrkStrawHits and the list of allowed states.
      TrkStrawHitStateVector(TSHV const& tshv,TSHSSV const& allowed);
      // copy constructor
      TrkStrawHitStateVector(TrkStrawHitStateVector const& other);
      TrkStrawHitStateVector& operator =(TrkStrawHitStateVector const& other);
      // set the state of the referenced TrkStrawHits.
      void setHitStates() const;
      // access the raw states
      TSHSV const & states() const { return _tshsv; }
      TSHSV & states() { return _tshsv; }
      // total # of states
      size_t nStates() const { return _nstates; }
      size_t nHits() const { return _tshsv.size(); }
      // increment or decrement the state, false return means not possible.  This is exhaustive
      bool increment();
      bool decrement();
      void reset();
      // return the sum of all penalties for the states of each TrkStrawHit
      static unsigned ipow(unsigned base, unsigned exp);
     private:
      TSHSV _tshsv; // vector of actua states
      TSHSSV _allowed; // allowed states
      std::vector<int> _state; // integer representation of the state
      unsigned _nstates; // cache of the number of states
      bool incrementState(); // increment state integers;
      bool decrementState(); // decrement state integers;
      void resetState(); // reset state integers;
      void intToVect();  // transform integer into state vector for each hit
      void vectToInt(); // transform state vector of each hit into integer
  };

}
#endif
