//
//  $Id: 
//  $Author: 
//  $Date: 
//
//  Original author Pavel Murat

#ifndef AlgorithmID_HH
#define AlgorithmID_HH

#include <utility>

namespace mu2e {

  class AlgorithmID {
  public:

    enum {
      TrkPatRecBit   = 0,
      CalPatRecBit   = 1
    };

    AlgorithmID(); 
    AlgorithmID(const AlgorithmID & p);

    ~AlgorithmID();

    AlgorithmID& operator = (const AlgorithmID & p);

    void clear() ;

    short   BestID () const { return _bestID ; }
    short   AlgMask() const { return _algMask; }

    void    Set       (short ID, short Mask) { _bestID = ID; _algMask = Mask; }

    void    SetBestID (short ID  ) { _bestID  = ID; }
    void    SetAlgMask(short Mask) { _algMask = Mask;}
    
  private:
    short   _bestID;			// bit ID identifing the best algorithm 
    short    _algMask;			// algorithm mask, for algorithms found this track
  };



} // end namespace mu2e


#endif
