#ifndef ITrackerGeom_ITLayerId_hh
#define ITrackerGeom_ITLayerId_hh

#include "ITrackerGeom/inc/SuperLayerId.hh"

namespace mu2e { 

struct ITLayerId{

        friend struct SuperLayerId;

public:

  ITLayerId():
    _sid(),
    _id(-1)
    {
  }
  
  ITLayerId( SuperLayerId *sid,
           int &id
           ):
    _sid(sid),
    _id(id)
  {
  }

  ~ITLayerId  (){
  }

  const SuperLayerId& getSuperLayerId() const {
    return *_sid;
  }

  const int& getSuperLayer() const {
    return _sid->_id;
  }

  const int getLayer() const{
    return _id;
  }

  bool operator==(const ITLayerId l) const{
    return ( *_sid == *(l._sid) && _id == l._id );
  }

  SuperLayerId *_sid;
  int _id;
  
};

inline std::ostream& operator<<(std::ostream& ost, 
                                const ITLayerId& l ){
  ost << "Layer Id: ("
      << l.getSuperLayerId() << " "
      << l._id
      << " )";
  return ost;
}

}
#endif /* ITrackerGeom_ITLayerId_hh */
