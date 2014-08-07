#ifndef MCDataProducts_CRVHit_hh
#define MCDataProducts_CRVHit_hh
//
// A CRV hit is a time stamp of a single PE, which was collected in one of the two fibers
// either on the positive or negative side of the counter;
//
// $Id: CRVHit.hh,v 1.1 2014/08/07 01:33:41 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//


namespace mu2e 
{
  class CRVHit
  {
    public:

    // Default c'tor needed by ROOT persistency.
    CRVHit(): _time(NAN), _fiberNumber(-1), _side(-1) {}

    CRVHit(double time, int fiberNumber, int side): _time(time), _fiberNumber(fiberNumber), _side(side) {}

    double time() const        {return _time;}
    int    fiberNumber() const {return _fiberNumber;}
    int    side() const        {return _side;}

    void   setTime(double time)            {_time=time;}
    void   setFiberNumber(int fiberNumber) {_fiberNumber=fiberNumber;}
    void   setSide(int side)               {_side=side;}

    //comparison operators needed for easy time ordering
    bool operator==( CRVHit const & rhs) const{return (_time == rhs._time);}
    bool operator< ( CRVHit const & rhs) const{return (_time <  rhs._time);}
    bool operator> ( CRVHit const & rhs) const{return (_time >  rhs._time);}

    private:

    double _time;
    int    _fiberNumber, _side;   //0 and 1 for both variables, 
                                  //where 0 is always the fiber or side which is located at the lower value 
                                  //at the x, y, or z axis 
  };
}

#endif /* MCDataProducts_CRVHit_hh */
