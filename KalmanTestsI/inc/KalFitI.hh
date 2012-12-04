//
// Object to perform BaBar Kalman fit
//
// $Id: KalFitI.hh,v 1.3 2012/12/04 00:51:26 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2012/12/04 00:51:26 $
//
#ifndef KalFitI_HH
#define KalFitI_HH

// framework
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTestsI/inc/TrkCellHit.hh"
// C++
#include <iostream>
#include <vector>

namespace mu2e 
{
//  class TrkStrawHit;

  class KalFitI : public KalFit
  {
  public:
// define the fit direction as downstream (towards positive Z) or upstream (towards negative Z).
    //enum fitDirection {downstream=0,upstream};
// define different t0 strategies.  Eventually t0 finding should be its own class
    //enum t0Strategy {external=0,median,histogram};
// parameter set should be passed in on construction
    explicit KalFitI(fhicl::ParameterSet const&);
    virtual ~KalFitI();
// main function: given a track definition, create a fit object from it
    virtual void makeTrack(KalFitResult& kdef);
    void makeExtrapolOutTrack(KalFitResult& kres);
    virtual const TrkVolume* trkVolume(trkDirection trkdir);

    // add a set of hits to an existing fit
    virtual void addHits(KalFitResult& kdef,std::vector<hitIndex> indices,bool active=true);
    void addHitsUnique(KalFitResult& kdef,std::vector<hitIndex> indices,bool active=true);

    void reActivateHitsbyTurn(KalFitResult& kdef);
    void reActivateHitsbyChi2(KalFitResult& kdef);

    void setMaxdist4Addhit(double maxd){_maxdist4Addhit=maxd;}

    bool useSingleCellDescr() { return _detailmat; }
    bool useDetailedWrSuppDescr() { return _detailWrSuppmat; }

// sim helper
    double _flt0;
    std::vector<double> _hitflt;

  private:
    // configuration parameters
    bool _doMinPrints;
    bool _detailmat;
    bool _detailWrSuppmat;
//    bool _ambigflip;
//    bool _updatet0;
//    double _t0tol;

    unsigned _maxfltdif;
    unsigned _maxMatReinter;
//    double _herr;

    double _maxdist4Addhit;

    //t0Strategy _t0strategy;
    bool _exacthitturn;

    double _momcuttoallowextrpl;
    // helper functions
    //virtual bool weedHits(KalFitResult& kdef);
    //virtual bool unweedHits(KalFitResult& kdef);
    virtual void fitTrack(KalFitResult& kdef);
    virtual void makeHits(KalFitResult& kdef,TrkT0 const& t0);

    virtual void makeMaterials(KalFitResult& kdef);
    void addHitMaterials(KalFitResult& kres, TrkCellHit* trkhit);

    bool fixHitTurn(KalFitResult& kres,TrkCellHit* trkhit);
    //virtual void initT0(TrkDef const& mytrk,TrkT0& t0);
    //virtual bool updateT0(KalFitResult& kdef);
    double findZFltlen(KalFitResult& kdef,double zval);
// general
//    static const double _vlight;
//    static const double _vdrift;
// 
  };
}
#endif
