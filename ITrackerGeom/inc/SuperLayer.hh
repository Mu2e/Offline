#ifndef SUPERLAYER_HH
#define SUPERLAYER_HH

#include <vector>
#include <deque>

#include "ITrackerGeom/inc/SuperLayerId.hh"
#include "ITrackerGeom/inc/SuperLayerInfo.hh"
#include "ITrackerGeom/inc/ITLayer.hh"


namespace mu2e {

class SuperLayer{

  friend class ITLayer;
  friend class ITracker;
  friend class ITrackerMaker;

public:

  // A free function, returning void, that takes a const Layer& as an argument.
  typedef void (*SuperLayerFunction)( const SuperLayer& s);

  SuperLayer();

  SuperLayer(SuperLayerId& id);

  SuperLayer(SuperLayerId&   id,
                  std::vector< boost::shared_ptr<ITLayer> > &layer);

  SuperLayer( int &id);

  SuperLayer( int &id,
                  std::vector< boost::shared_ptr<ITLayer> > &layer);

  ~SuperLayer ();
 
  const SuperLayerId& Id() const { return _id;}
  
  int nLayers() const { return _nLayers; }

  boost::shared_ptr<ITLayer> getLayer( int n ) const throw(cet::exception) {
    if (n>=0 && n<_nLayers) return _layers.at(n);
    else throw  cet::exception("GEOM")<< "Layer number: "<< n <<" not present in "<<_id;
  }
  
  boost::shared_ptr<ITLayer> getLayer( ITLayerId& id ) const {
    return getLayer(id._id);
  }
  
  const std::vector< boost::shared_ptr<ITLayer> >& getLayers() const {
    return _layers;
  }

//  SuperLayer& operator=(const SuperLayer &sl) {
//          if (this!=&sl) {
//                  _id = sl.Id();
//                  _nLayers = sl.nLayers();
//                  _layers = sl.
//          }
//  }

protected:

  SuperLayerId _id;

  int _nLayers;

  std::vector< boost::shared_ptr<ITLayer> > _layers;

  void addLayer(ITLayer *itl){
          _layers.push_back(boost::shared_ptr<ITLayer>(itl));
          _nLayers++;
  }

};
}

#endif /*SUPERLAYER_HH*/
