#ifndef TrackerGeom_ManifoldId_hh
#define TrackerGeom_ManifoldId_hh

#include <ostream>

#include "DataProducts/inc/PanelId.hh"

namespace mu2e {

  struct ManifoldId{

  public:

    ManifoldId():
      _pnlid(),
      _manifold(-1){
    }

    ManifoldId(PanelId panel,
               int manifold
               ):
      _pnlid(panel),
      _manifold(manifold)
    {}

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    PanelId getPanelId() const {
      return _pnlid;
    }

    int getManifold() const {
      return _manifold;
    }

    bool operator== (const ManifoldId s) const{
      return (_pnlid == s._pnlid && _manifold == s._manifold);
    }

  private:

    // Member variables
    PanelId _pnlid;
    int _manifold;

  };

  inline std::ostream& operator<<(std::ostream& ost,
                                  const ManifoldId& m ){
    ost << m.getPanelId() << " " << m.getManifold();
    return ost;
  }

} // end namespace mu2e

#endif /* TrackerGeom_ManifoldId_hh */
