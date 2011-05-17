#ifndef ITrackerGeom_Wall_hh
#define ITrackerGeom_Wall_hh

#include <vector>
#include <boost/shared_ptr.hpp>

//#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/Transform3D.h"
//#include "ITrackerGeom/inc/ITLayer.hh"

namespace mu2e {

class Wall{

        friend class ITrackerMaker;

public:

        // A free function, returning void, that takes a const Layer& as an argument.
        typedef void (*WallFunction)( const Wall& s);

        enum Walltype {undefined=-1, inner, outer, endcap};

        Wall();

        Wall(Walltype wt=Wall::undefined);

        Wall(const Wall &wl);

        ~Wall (){}

        //        CLHEP::Hep3Vector getC()                const { return _c; }

        HepGeom::Transform3D getPos() const { return _pos; }

        Walltype getType()            const { return _type; }

        int getNShells()              const { return _nShells; }

        float getDPhi()               const { return _pDPhi; }

        float getDTheta()             const { return _pDTheta; }

        float getDz()                 const { return _pDz; }

        std::string getName()         const { return _name; }

        float getRmax()               const { return _pRmax; }

        float getRmin()               const { return _pRmin; }

        float getSPhi()               const { return _pSPhi; }

        float getSTheta()             const { return _pSTheta; }

        float getTotalThickness()     const {
                return _totalThickness;
        }

        boost::shared_ptr<std::vector<std::string> > getMaterialsName() const {
                return _materialsName;
        }

        boost::shared_ptr<std::vector<double> > getThicknesses() const {
                return _thicknesses;
        }

        void addMaterials(int &wShellNumber, std::vector<std::string> *wShellsMatName, std::vector<double> *wShellsThicknesses) throw(cet::exception);
        /*{
                _nShells=wShellNumber;
                if ( _nShells!=wShellsMatName.size() && _nShells!=wShellsThicknesses.size() )
                        throw cet::exception("GEOM")<< "Error in Configuration file! There is a disagreement between the vectors dimensions of a ITracker wall.\n";

                _materialsName.reset(&wShellsMatName);
                _thicknesses.reset(&wShellsThicknesses);
                _totalThickness = 0.0;
                for (int is = 0; is < _nShells; ++is) {
                        _totalThickness += wShellsThicknesses.at(is);
                }
        }*/

        Wall& operator=(const Wall &wl);
        /*{
                if (this!=&wl) {
                        _type=wl.getType();
                        _nShells=wl.getShells();
                        _totalThickness=wl.getTotalThickness();
                        _pRmin=wl.getRmin();
                        _pRmax=wl.getRmax();
                        _pSPhi=wl.getSPhi();
                        _pDPhi=wl.getDPhi();
                        _pSTheta=wl.getSTheta();
                        _pDTheta=wl.getDTheta();
                        _pDz=wl.getDz();
                        _name=wl.getName();
                        _c=wl.getC();
                        _materialsName=wl.getMaterialsName();
                        _thicknesses=wl.getThicknesses();
                }
        }*/

protected:

        Walltype     _type;

        float        _pRmin;
        float        _pRmax;
        float        _pSPhi;
        float        _pDPhi;
        float        _pSTheta;
        float        _pDTheta;
        float        _pDz;
        std::string  _name;
        //CLHEP::Hep3Vector  _c;
        HepGeom::Transform3D _pos;

private:
        int                                          _nShells;
        float                                        _totalThickness;
        boost::shared_ptr<std::vector<std::string> > _materialsName;
        boost::shared_ptr<std::vector<double> >      _thicknesses;

};

} //namespace mu2e

#endif /* ITrackerGeom_Wall_hh */
