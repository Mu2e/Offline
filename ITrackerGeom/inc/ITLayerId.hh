#ifndef ITrackerGeom_ITLayerId_hh
#define ITrackerGeom_ITLayerId_hh

#include "ITrackerGeom/inc/SuperLayerId.hh"

namespace mu2e {

class ITLayerId;
inline std::ostream& operator<<(std::ostream& ost,
                                const ITLayerId& l );

class ITLayerId{

        friend struct SuperLayerId;

public:

  ITLayerId():
    _sid(),
    _id(-1)
  {
  }

  ITLayerId( SuperLayerId *sid,
             int id
             ):
    _sid(sid),
    _id(id)
  {
  }

  // Use compiler-generated copy c'tor, copy assignment, and d'tor

  const SuperLayerId& getSuperLayerId() const {
    return *_sid;
  }

  int getSuperLayer() const {
    return _sid->getSuperLayer();
  }

  int getLayer() const{
    return _id;
  }

  bool operator==(const ITLayerId& l) const{
    return ( *_sid == *(l._sid) && _id == l._id );
  }

  friend std::ostream& operator<<(std::ostream& ost,
                                  const ITLayerId& l ){
    ost << "Layer Id: ("
        << l.getSuperLayerId() << " "
        << l._id
        << " )";
    return ost;
  }

private:

  SuperLayerId *_sid;
  int _id;

};

}
#endif /* ITrackerGeom_ITLayerId_hh */
