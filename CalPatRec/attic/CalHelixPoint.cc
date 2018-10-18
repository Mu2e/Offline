//

#include "CalPatRec/inc/CalHelixPoint.hh"


namespace mu2e {

//-----------------------------------------------------------------------------
// CalHelixPoint
//-----------------------------------------------------------------------------
  // statics
  //  double       CalHelixPoint::_efac(1.0);
  StrawHitFlag CalHelixPoint::_useflag;

  CalHelixPoint::CalHelixPoint(size_t                  index, 
			       ComboHit const&         sh, 
			       // StrawHitPosition const& shp, 
			       // Straw const&            straw,
			       StrawHitFlag const&     flag) :
    // _ind (index), 
    // _pos (shp.posCLHEP()), 
    // _flag(flag),
    // _wdir(straw.getDirection()),
    // _straw(&straw),
    // _strawhit(&sh),
    // _perr(_efac*shp.posRes(StrawHitPosition::wire)),
    // _rerr(_efac*shp.posRes(StrawHitPosition::trans))
    // _perr(shp.posRes(StrawHitPosition::wire)),
    // _rerr(shp.posRes(StrawHitPosition::trans)),
    ComboHit(sh),
    _dzFromSeed(0),
    _drFromPred(0),
    _xyWeight(0),
    _zphiWeight(0)    
  {
    _ncombo  = 1;
    _pind[0] = index;
    _flag    = flag;

    static const XYZVec _zdir(0.0,0.0,1.0);
    _sdir = _zdir.Cross(_wdir);
  }

  // CalHelixPoint::CalHelixPoint(size_t                   ind, 
  // 			       CLHEP::Hep3Vector const& pos, 
  // 			       CLHEP::Hep3Vector const& wdir, 
  // 			       double                   werr, 
  // 			       double                   serr) :
  //   _ind(ind),
  //   _pos(pos),
  //   //    _phi(_pos.phi()),
  //   _flag(),
  //   _wdir(wdir),
  //   _sdir(wdir.y(),-wdir.x(),0.0),
  //   _straw(0),
  //   _strawhit(0),
  //   // _perr(_efac*werr),
  //   // _rerr(_efac*serr)
  //   _perr(werr),
  //   _rerr(serr),
  //   _dzFromSeed(0),
  //   _drFromPred(0),
  //   _testDzFromSeed(0),
  //   _testDrFromPred(0),
  //   _xyWeight(0),
  //   _zphiWeight(0)    
  // {
  // }
  
  CalHelixPoint::CalHelixPoint(const CalHelixPoint& Copy){
    *this    = Copy;
    _useflag = Copy._useflag;	
    // _ind     = Copy._ind;		
    // _pos     = Copy._pos;		
    // _phi     = Copy._phi;	        
    // _flag    = Copy._flag;		
    // _wdir    = Copy._wdir;		
    _sdir    = Copy._sdir;           
    			
    			
    // _straw    = Copy._straw;          
    // _strawhit = Copy._strawhit;       
    // _perr     = Copy._perr;
    // _rerr     = Copy._rerr;
    
    _dzFromSeed = Copy._dzFromSeed;     
    _drFromPred = Copy._drFromPred;     

    _xyWeight       = Copy._xyWeight;
    _zphiWeight     = Copy._zphiWeight;    
  }

  bool CalHelixPoint::use() const { 
    return !_flag.hasAnyProperty(_useflag);
  }

  bool CalHelixPoint::stereo() const {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return _flag.hasAllProperties(stereo);
  }

  // void CalHelixPoint::setUse(bool use) {
  //   static StrawHitFlag other(StrawHitFlag::other);
  //   if(!use)
  //     _flag.merge(other);
  //   else
  //     _flag.clear(other);
  // }

  void CalHelixPoint::setOutlier(){
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    _flag.merge(outlier);
  }

  bool CalHelixPoint::isOutlier() const {
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    return _flag.hasAllProperties(outlier);
  }

  bool CalHelixPoint::isCalosel() const {
    static StrawHitFlag calosel(StrawHitFlag::calosel);
    return _flag.hasAllProperties(calosel);
  }

}
