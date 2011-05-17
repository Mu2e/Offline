#ifndef ITRACKER_HH
#define ITRACKER_HH

#include <deque>
#include <vector>
#include <map>
#include <iterator>

#include <boost/shared_array.hpp>

#include "ITrackerGeom/inc/SuperLayer.hh"
#include "ITrackerGeom/inc/Wall.hh"
//#include "GeometryService/inc/Detector.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "ITrackerGeom/inc/CellGeometryHandle.hh"

namespace mu2e {

class ITracker: public Tracker {

        friend class ITrackerMaker;

public:
        ITracker();
        ~ITracker() {}

        virtual std::string name() const { return "ITracker";}


        enum GeomType          { Hexagonal=2, Square };
        enum EnCapType         { Plane, Spherical };

        double r0()            const { return _r0;}
        double z0()            const { return _z0;}
        double rOut()          const { return _rOut;}
        int nSWire()           const { return _nSWire;}
        int nSDeltaWire()      const { return _nSDeltaWire;}
        int nRing()            const { return _nRing;}

        std::string extFile()  const { return _extFile; }
        bool isExternal()      const { return _isExternal; }

        int nSuperLayers()     const { return _nSuperLayers; }

        double zHalfLength()   const { return _zHalfLength;}

        double maxEndCapDim()  const { return _max_EndCap_dim; }

        GeomType geomType()    const { return _geomType; }

        EnCapType endcapType() const { return _endcapType; }

        int getNWalls()        const { return _nWalls; }

        bool displayGasLayer() const { return _displayGasLayer; }
        bool displayWires()    const { return _displayWires; }

        CellGeometryHandle* getCellGeometryHandle() const { return _cellhnd.get(); }

        SuperLayer* getSuperLayer(int n) const throw(cet::exception);

        boost::shared_array<SuperLayer> getSuperLayersArray() const {
                return _sprlr;
        }

        const Straw& getStraw ( const StrawId& sid )      const throw(cet::exception) {
                throw cet::exception("GEOM")<< "Fake method \"getStraw ( StrawId )\", not used for the ITracker";
                return fakeStraw;
        }
//        const Straw& getStraw ( StrawIndex i )      const { return fakeStraw; }
        const Straw&  getStraw ( StrawIndex i )            const {
                getCellGeometryHandle()->SelectCellDet(i.asUint());
                return dynamic_cast<Straw&>( *(getCellGeometryHandle()->GetITCell().get()) );
        }
        const std::deque<Straw>& getAllStraws()           const throw(cet::exception) {
                throw cet::exception("GEOM")<< "Fake method \"getAllStraws()\", not used for the ITracker";
                return fakeStrawDeq;
        }
        const std::vector<StrawDetail>& getStrawDetails() const throw(cet::exception) {
                throw cet::exception("GEOM")<< "Fake method \"getStrawDetails()\", not used for the ITracker";
                return fakeStrawVec;
        }

//        const boost::shared_ptr<Wall> getWall(int n) throw(cet::exception);

        boost::shared_ptr<std::multimap<Wall::Walltype,boost::shared_ptr<Wall> > > getWalls() const {
                return _walls;
        }

protected:

        // Nominal values.
        // _r0 = Nominal radius of the center of the sector.
        // _z0 = position of the center of the tracker relative to the origin
        //       of the Mu2e coordinate system.
        double _r0;
        double _z0;
        int _nSWire;
        int _nSDeltaWire;
        int _nRing;

        // Outer radius of a logical volume that will just contain the entire tracker.
        double _rOut;

        // Name of external gdml geometry file description.
        std::string _extFile;
        bool _isExternal;

        int _nSuperLayers;
        int _nWalls;

        double _zHalfLength;
        double _max_EndCap_dim;

        //Cell geometry type: 2:Hexagonal, 3:Square
        GeomType _geomType;

        //EndCap shape type: 0 plane, 1 spherical
        EnCapType _endcapType;

        //Allow to display the gas (and/or every wires inside gas) inside the chamber
        bool _displayGasLayer;
        bool _displayWires;

        boost::shared_array<SuperLayer> _sprlr;

        boost::shared_ptr<std::multimap<Wall::Walltype,boost::shared_ptr<Wall> > > _walls;

        void addWall(Wall *wall);

        std::auto_ptr<CellGeometryHandle> _cellhnd;

private:
        std::multimap <Wall::Walltype,boost::shared_ptr<Wall> >::iterator _walls_it;
        int _lastSeenWall; //last extracted wall from the _walls container using the getWall method

        const Straw fakeStraw;
        const std::deque<Straw> fakeStrawDeq;
        const std::vector<StrawDetail>  fakeStrawVec;
};

} //namespace mu2e

#endif /*ITRACKER_HH*/
