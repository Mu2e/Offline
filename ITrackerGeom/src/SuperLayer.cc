#include "ITrackerGeom/inc/SuperLayer.hh"

#ifndef __CINT__ 

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e {

SuperLayer::SuperLayer():
    		_id(SuperLayerId()),
    		_nLayers(0),
    		_layers(*(new std::vector< boost::shared_ptr<ITLayer> >(0)))
    		{
    		}

SuperLayer::SuperLayer(SuperLayerId& id):
			_id(id),
			_nLayers(0),
			_layers(*(new std::vector< boost::shared_ptr<ITLayer> >(0)))
			{
			}

SuperLayer::SuperLayer(SuperLayerId& id, std::vector< boost::shared_ptr<ITLayer> > &layers):
			_id(id),
			_layers(layers)
			{
	_nLayers = layers.size();
			}

SuperLayer::SuperLayer( int& id):
			_id(SuperLayerId(id)),
			_nLayers(0),
			_layers(*(new std::vector< boost::shared_ptr<ITLayer> >(0)))
			{
			}

SuperLayer::SuperLayer( int& id, std::vector< boost::shared_ptr<ITLayer> > &layers):
			_id(SuperLayerId(id)),
			_layers(layers)
			{
	_nLayers = layers.size();
			}

SuperLayer::~SuperLayer(){
	//	 for ( std::vector<ITLayer*>::iterator j=_layers.begin(); j != _layers.end(); j++){
	//		 delete (*j);
	//	 }
}

} // namespace mu2e 

#endif

