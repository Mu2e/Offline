#include "ITrackerGeom/inc/ITracker.hh"

using namespace std;

namespace mu2e {

ITracker::ITracker() {
        _r0              =0.0;
        _z0              =0.0;
        _nSWire          =0;
        _nSDeltaWire     =0;
        _nRing           =0;
        _isDumbbell      =false;
        _rOut            =0.0;
        _extFile         ="";
        _isExternal      =false;
        _nSuperLayers    =0;
        _nWalls          =0;
        _zHalfLength     =0.0;
        _max_EndCap_dim  =0.0;
        _geomType        =ITracker::Hexagonal;
        _endcapType      =ITracker::Plane;
        _displayGasLayer =false;
        _displayWires    =false;
        _zZonesLimits.reset(new double[2]);
        _zZonesLimits.get()[0] =0.0;
        _zZonesLimits.get()[1] =0.0;
}

SuperLayer* ITracker::getSuperLayer(int n) const throw(cet::exception) {
        if (n>=0 && n< _nSuperLayers){
                return &(_sprlr[n]);
        }
        else throw cet::exception("GEOM")<< "Super Layer number: "<< n <<" not present";
}

//        const boost::shared_ptr<Wall> ITracker::getWall(int n) throw(cet::exception) {
//                if (n>=0 && n< _nWalls){
//                        //_walls_it --;//+= (n-_lastSeenWall);
//                        advance (_walls_it,n-_lastSeenWall);
//                        _lastSeenWall = n;
//                        return _walls_it->second;
//                }
//                else throw cet::exception("GEOM")<< "Wall number: "<< n <<" not present";
//        }

void ITracker::addWall(Wall *wall){
        if (_nWalls==0) _walls.reset(new std::multimap<Wall::Walltype,boost::shared_ptr<Wall> >() );
        _walls->insert(std::pair<Wall::Walltype,boost::shared_ptr<Wall> >(wall->getType(),boost::shared_ptr<Wall>(wall)) );
        _nWalls++;
        if (_nWalls==1){
                _walls_it=_walls->begin();
                _lastSeenWall=0;
        }
}

} // namespace mu2e
