//
// Fast Patter recognition Data type and method of general use
//
// $Id: FastPatRecoUtilsAndDataDef.hh,v 1.4 2012/05/22 06:37:04 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/05/22 06:37:04 $
//
// Original author G. Tassielli
//
#ifndef FastPatRecoUtilsAndDataDef_HH
#define FastPatRecoUtilsAndDataDef_HH  

// C++ includes.
#include <set>
//#include <unordered_set>
#include <map>
#include <algorithm>
#include <limits>
#include <cmath>

#include <boost/shared_ptr.hpp>

// Framework includes.
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"

#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "FastPatternReco/inc/KarimakiCircle.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/TrackerHitByID.hh"

using namespace std;

namespace mu2e {

typedef std::set<size_t> isHitIDUsed;

//----------------- for TTtracker -----------------

typedef std::multimap<unsigned int, size_t, less<unsigned int> > stbrel;                                        //"Straw-histogram bin" relation (navigation)

typedef std::multimap<unsigned short, unsigned short, less<unsigned short> > ptcsctrrel;                        //"Station Pitch-Station" relation
typedef std::map<unsigned short, ptcsctrrel> ptcbnrel;                                                          //container for each row "Station Pitch-Station" relation
typedef std::set<unsigned short, less<unsigned short> > goodsttnrel;                                                               //"good Station" that has at least one vote in both scan direction per row
typedef std::map<unsigned short, goodsttnrel> goodsctrstatnrel;                                                 //container for each row "good Station"
typedef std::multimap<unsigned short, std::pair<size_t, unsigned short>, less<unsigned short> > ptcclsbrel;     //"Station Pitch-row Cluster" relation ("row Cluster" = cluster along each row in the Station-Sector map)
typedef std::multimap<unsigned short, ptcclsbrel, less<unsigned short> > ptcmaprowclsrel;                       //container of the "Station Pitch-row Cluster" relation for each sector (row should be see like an alias for sector)
typedef std::pair<size_t, size_t> rwclclcpl;                                                                    //"row Cluster-row Cluster" coupling along a row
typedef std::map<rwclclcpl, std::vector<unsigned short> > avptcclscpl;                                          //list of all available station pitch for each "row Cluster-row Cluster" pair
typedef std::multimap<unsigned short, avptcclscpl , less<unsigned short> > rwavptcclscpl;                       //container of the data for "row Cluster-row Cluster" pair for each row


typedef art::Ptr<StrawHit> StrawHitPtr;
typedef std::multimap<unsigned int, StrawHitPtr, less<unsigned int> > stMaprel;
typedef art::Ptr<TrackerHitTimeCluster> TrackerHitTimeClusterPtr;


//----------------- for ITtracker -----------------

int nCellPerLayer = 0;

typedef std::vector< std::pair<int, size_t> > hitsInClsID; //first = cell ID in Layer, second = hit ID in the event

struct closHitClust{

        closHitClust():
                _nHit(0),
                _maxCellID(0),
                _minCellID(0),
                _centerCellID(0)
        {}

        closHitClust(unsigned short cellID, size_t hitIdx):
                _nHit(1),
                _maxCellID(cellID),
                _minCellID(cellID),
                _centerCellID(cellID)
        {
                _hitIdx.push_back(std::pair<int, size_t>(cellID,hitIdx) );
        }

        //int maxCellID() { return (_hitIdx.back().first); }
        //int minCellID() { return (_hitIdx.front().first); }
        unsigned int _nHit;
        int _maxCellID;
        int _minCellID;
        int          _centerCellID;
        hitsInClsID _hitIdx; //first = cell ID in Layer, second = hit ID in the event

        void addHit(unsigned short cellID, int distStartingHit, size_t hitIdx) {
                ++_nHit;
                if (distStartingHit>0) {
                        _hitIdx.push_back(std::pair<int, size_t>(cellID,hitIdx));
                } else {
                        _hitIdx.insert(_hitIdx.begin(),std::pair<int, size_t>(cellID,hitIdx));
                }
                _maxCellID=_hitIdx.back().first;
                _minCellID=_hitIdx.front().first;
                _centerCellID=_hitIdx.back().first-_nHit/2;
                if (_nHit%2==0) ++_centerCellID;
                if (_centerCellID<0) _centerCellID+=nCellPerLayer;

        }

        bool operator==( closHitClust const& comp) const {
                return ( _nHit==comp._nHit &&
                                _centerCellID==comp._centerCellID &&
                                _hitIdx.front().second==comp._hitIdx.front().second &&
                                _hitIdx.back().second==comp._hitIdx.back().second );
        }

        bool operator!=( closHitClust const& comp) const {
                return ( _nHit!=comp._nHit ||
                                _centerCellID!=comp._centerCellID ||
                                _hitIdx.front().second!=comp._hitIdx.front().second ||
                                _hitIdx.back().second!=comp._hitIdx.back().second );
        }

        bool operator<( closHitClust const& comp) const {
                return ( _nHit<comp._nHit );
        }

        bool operator>( closHitClust const& comp) const {
                return ( _nHit>comp._nHit );
        }

};

bool operator >( const boost::shared_ptr<mu2e::closHitClust> &c1, const boost::shared_ptr<mu2e::closHitClust> &c2 )
{
        if ( c1->_nHit > c2->_nHit ) { return true; }
        if ( ( c1->_nHit == c2->_nHit ) &&  ( c1->_centerCellID > c2->_centerCellID ) ) { return true; }
        return false;
        //return ( ( c1->_nHit >= c2->_nHit ) && ( c1->_centerCellID > c2->_centerCellID ) );
}

/*bool operator <( const boost::shared_ptr<mu2e::closHitClust> &c1, const boost::shared_ptr<mu2e::closHitClust> &c2 )
  {
          return ( ( c1->_nHit <= c2->_nHit ) && ( c1->_centerCellID < c2->_centerCellID ) );
  }*/

typedef std::set< boost::shared_ptr<mu2e::closHitClust>, std::greater< boost::shared_ptr<mu2e::closHitClust> > > closClstcol;
typedef std::map<int, closClstcol, std::greater<int> > closClinRadLay;   //"closest hit Clusters" for each Radial Layer (key=Radial Layer Id)

struct ptrHitInClosCust {

        ptrHitInClosCust() {}

        ptrHitInClosCust(hitsInClsID::const_iterator &hitsInClsID_it, closClstcol::iterator &closClstcol_it, closClinRadLay::iterator &closClinRadLay_it) {
                _hitsInClsID_it    = hitsInClsID_it;
                _closClstcol_it    = closClstcol_it;
                _closClinRadLay_it = closClinRadLay_it;
        }

        void setHitInClosCust(hitsInClsID::const_iterator &hitsInClsID_it, closClstcol::iterator &closClstcol_it, closClinRadLay::iterator &closClinRadLay_it) {
                _hitsInClsID_it    = hitsInClsID_it;
                _closClstcol_it    = closClstcol_it;
                _closClinRadLay_it = closClinRadLay_it;
        }

        int const & getRadLayID () const {
                return _closClinRadLay_it->first;
        }

        mu2e::closHitClust const & getHitClust () const {
                return *(*_closClstcol_it);
        }

        int const & getinLayerCellID () const {
                return _hitsInClsID_it->first;
        }

        size_t const & getinEventHitID () const {
                return _hitsInClsID_it->second;
        }

        bool operator==( ptrHitInClosCust const& comp) const {
                return ( _hitsInClsID_it->second==comp._hitsInClsID_it->second &&
                                _hitsInClsID_it->first==comp._hitsInClsID_it->first &&
                                _closClinRadLay_it->first==comp._closClinRadLay_it->first );
        }

        bool operator<( ptrHitInClosCust const& comp) const {
                return ( _hitsInClsID_it->second<comp._hitsInClsID_it->second );
        }

        friend std::ostream& operator<< ( std::ostream& ost,
                                  const ptrHitInClosCust& gtd_ ){
                ost<<"Cell "<<gtd_._hitsInClsID_it->first<<" RadLayer "<<gtd_._closClinRadLay_it->first<<" hit n. "<<gtd_._hitsInClsID_it->second<<
                                " cluster size "<<(*gtd_._closClstcol_it)->_nHit<<" minCell "<<(*gtd_._closClstcol_it)->_minCellID<<" maxCell "<<
                                (*gtd_._closClstcol_it)->_maxCellID<<" centralCell "<<(*gtd_._closClstcol_it)->_centerCellID<<endl;

                return ost;
        }

private :
        hitsInClsID::const_iterator _hitsInClsID_it;
        closClstcol::iterator _closClstcol_it;
        closClinRadLay::iterator _closClinRadLay_it;

};

struct confMapPoint{

        confMapPoint():
                _u(0.0),
                _v(0.0),
                _rSq(0.0)
        {}

        confMapPoint(double u, double v, double invRSq, double sigmaC):
                _u(u),
                _v(v)
        //_rSq(rSq)
        {
                _errU = sigmaC*sqrt( 2.0*invRSq*invRSq + 8.0*u*u );
                _errV = sigmaC*sqrt( 2.0*invRSq*invRSq + 8.0*v*v );
                _rSq = 1.0/invRSq;
        }

        double _u;
        double _v;
        double _rSq;
        double _errU;
        double _errV;

        bool operator==( confMapPoint const& comp) const {
                return ( _u==comp._u && _v==comp._v && _rSq==comp._rSq );
        }

};

//struct confMapDraw{
//
//        confMapDraw():
//                _cmap(0x0),
//                _cmapCHT(0x0)
//        {}
//
//        confMapDraw( TGraphErrors *cmap, TH2F *cmapCHT):
//                _cmap(cmap),
//                _cmapCHT(cmapCHT)
//        {}
//
//        ~confMapDraw() {
//                //cout<<" pointers "<<_cmap<<" "<<_cmapCHT<<endl;
//                if (_cmap!=0x0) delete _cmap;
//                if (_cmapCHT!=0x0) delete _cmapCHT;
//        }
//
//        //TH2F *_cmap;
//        TGraphErrors *_cmap;
//        TH2F *_cmapCHT;
//
//};

struct CHTVotArr {

        CHTVotArr(unsigned int nBinX=2513, unsigned int nBinY=500, float minX=-CLHEP::pi, float maxX=CLHEP::pi, float minY=-0.25, float maxY=0.25, unsigned int thr=10):
                _nBinX(nBinX),
                _nBinY(nBinY),
                _minX(minX),
                _maxX(maxX),
                _minY(minY),
                _maxY(maxY),
                _THR(thr)
        {
                _binWdtX  = (_maxX-_minX)/((float)_nBinX);
                _binWdtY  = (_maxY-_minY)/((float)_nBinY);
                _ovrfluID = _nBinX*_nBinY;
                _dataArr = new unsigned int*[_nBinX];
                for (unsigned int i=0; i<_nBinX; ++i) {
                        _dataArr[i] = new unsigned int[_nBinY];
                        for (unsigned int j=0; j<_nBinY; ++j) {
                                _dataArr[i][j]=0;
                        }
                }
        }

        ~CHTVotArr() {
                for (int i=0; i<_nBinX; ++i) {
                        delete [] _dataArr[i];
                }
                delete [] _dataArr;
        }

        void setTHR(unsigned int thr) {
                _THR = thr;
                if (_overTHRs.size()>0) {
                        _overTHRs.clear();
                        for (unsigned int i=0; i<_nBinX; ++i) {
                                for (unsigned int j=0; j<_nBinY; ++j) {
                                        if (_dataArr[i][j]>_THR){
                                                unsigned int tmpuID=_ovrfluID;
                                                xyTouID(i, j, tmpuID);
                                                _overTHRs.insert( tmpuID );
                                        }
                                }
                        }
                }
        }

        bool xyTouID(double x, double y, unsigned int &outuID) {
                unsigned int tmpX = (unsigned int)((x-_minX)/_binWdtX);
                unsigned int tmpY = (unsigned int)((y-_minY)/_binWdtY);
                outuID=_ovrfluID;
                return xyTouID (tmpX, tmpY, outuID);
        }

        bool xyTouID(unsigned int tmpX, unsigned int tmpY, unsigned int &outuID) {
                if ( tmpX<_nBinX && tmpY<_nBinY ) {
                        outuID = (tmpX*_nBinY + tmpY);
                        return true;
                }
                else { return false; }
        }

        bool uIDToxy(unsigned int uID, unsigned int &outX, unsigned int &outY) {
                if ( uID<_ovrfluID ) {
                        outY = uID%_nBinY;
                        outX = (uID-outY)/_nBinY;
                        return true;
                }
                else { return false; }
        }

        unsigned int Fill(float x, float y){
                unsigned int tmpX = (unsigned int)((x-_minX)/_binWdtX);
                unsigned int tmpY = (unsigned int)((y-_minY)/_binWdtY);
                unsigned int tmpuID=_ovrfluID;
                if ( xyTouID(tmpX, tmpY, tmpuID) ) {
                        ++_dataArr[tmpX][tmpY];
                        if (_dataArr[tmpX][tmpY]>_THR){
                                _overTHRs.insert( tmpuID );
                        }
                }
                return tmpuID;
        }

        unsigned int get(unsigned int tmpX, unsigned int tmpY) {
                if ( tmpX<_nBinX && tmpY<_nBinY ) {
                        return _dataArr[tmpX][tmpY];
                }
                else { return 0; }
        }

        unsigned int get(unsigned int uID) {
                unsigned int tmpX=0, tmpY=0;
                if ( uIDToxy( uID, tmpX, tmpY) ) {
                        return _dataArr[tmpX][tmpY];
                }
                else { return 0; }
        }

        std::set< unsigned int > _overTHRs; //list of the unique ID for the points over threshold

        float & getBinWdtX() { return _binWdtX; }
        float & getBinWdtY() { return _binWdtY; }
        float & getMinX()    { return _minX; }
        float & getMinY()    { return _minY; }

private:
        float _binWdtX;
        float _binWdtY;
        unsigned int _nBinX;
        unsigned int _nBinY;
        float _minX;
        float _maxX;
        float _minY;
        float _maxY;
        unsigned int _THR;
        unsigned int **_dataArr;

        unsigned int _ovrfluID;
};

typedef std::multimap<unsigned int , ptrHitInClosCust > CHTVotArr_HitPtrrel;

struct ClosClustCHTVot {

        ClosClustCHTVot():
                _meanX(0.0),
                _meanY(0.0),
                _sigmaX(0.0),
                _sigmaY(0.0),
                nHit(0.0),
                tmpData(0.0),
                tmpMeanX(0.0),
                tmpMSX(0.0),
                tmpMeanY(0.0),
                tmpMSY(0.0)
        {}

        ~ClosClustCHTVot() {}

        void addhit(unsigned int tmpX, unsigned int tmpY, unsigned int tmpuID, unsigned int iHitMult ) {
                _listCHTVotIDs.insert(tmpuID);
                nHit     += (float)iHitMult;
                tmpData   = (float)tmpX + 0.50000;
                tmpMeanX += ((float)iHitMult)*tmpData;
                tmpMSX   += ((float)iHitMult)*tmpData*tmpData;
                _meanX    = tmpMeanX/nHit;
                invSqrtNhit_1 = 1.0/sqrt( nHit-1.0 );
                tmpBuff   = tmpMSX/nHit - _meanX*_meanX;
                if (tmpBuff>0.0) {
                        _sigmaX = sqrt( tmpBuff * nHit ) * invSqrtNhit_1;
                }else {
                        _sigmaX = 0.288675135 * invSqrtNhit_1;   // 1/sqrt(12) * 1/sqrt(n-1)
                }
                tmpData   = (float)tmpY + 0.50000;
                tmpMeanY += ((float)iHitMult)*tmpData;
                tmpMSY   += ((float)iHitMult)*tmpData*tmpData;
                _meanY    = tmpMeanY/nHit;
                tmpBuff   = tmpMSY/nHit - _meanY*_meanY;
                if (tmpBuff>0.0) {
                        _sigmaY = sqrt( tmpBuff * nHit ) * invSqrtNhit_1;
                } else {
                        _sigmaY = 0.288675135 * invSqrtNhit_1;   // 1/sqrt(12) * 1/sqrt(n-1)
                }
        }

        ClosClustCHTVot & operator= ( ClosClustCHTVot const& tmpClus ) {
                _meanX         = tmpClus._meanX;
                _meanY         = tmpClus._meanY;
                _sigmaX        = tmpClus._sigmaX;
                _sigmaY        = tmpClus._sigmaY;
                nHit           = tmpClus.nHit;
                tmpData        = tmpClus.tmpData;
                tmpMeanX       = tmpClus.tmpMeanX;
                tmpMSX         = tmpClus.tmpMSX;
                tmpMeanY       = tmpClus.tmpMSY;
                tmpMSY         = tmpClus.tmpMSY;
                _listCHTVotIDs = tmpClus._listCHTVotIDs;
        }

        std::set< unsigned int > _listCHTVotIDs; //list of the unique ID of the points in the voting array

        float _meanX;
        float _meanY;
        float _sigmaX;
        float _sigmaY;

private:
        float nHit;
        float tmpData;
        float tmpMeanX;
        float tmpMSX;
        float tmpMeanY;
        float tmpMSY;

        float tmpBuff;
        float invSqrtNhit_1;
};

struct SimpleCircle2D {

        SimpleCircle2D():
                _radius(0.0),
                _sigmaRad(0)
        {
                _center.set(0.0,0.0);
                _sigmaCenter[0] = _sigmaCenter[1] = 0.0;
        }

        ~SimpleCircle2D() {}

        void SetCenter (float x, float y, float sigmaX=0.0, float sigmaY=0.0) {
                _center.set(x,y);
                _sigmaCenter[0] = sigmaX;
                _sigmaCenter[1] = sigmaY;
        }

        void SetRadius (float rad, float sigmaRad=0.0) {
                _radius = rad;
                _sigmaRad = sigmaRad;
        }

        bool isCompatibleByCmapWith( SimpleCircle2D const& tmpCirc, float sigmaLevel = 1.0 ) {
                float tmpSigma = sqrt(_sigmaRad*_sigmaRad + tmpCirc._sigmaRad*tmpCirc._sigmaRad);
                float tmpDist  = _radius - tmpCirc._radius;
                tmpSigma *= sigmaLevel;
                if ( tmpDist>-tmpSigma && tmpDist<tmpSigma) {
                        //cout<<"Compatible by R"<<endl;
                        float dx2, dy2;
                        dx2 = _center.x()-tmpCirc._center.x();
                        dx2*=dx2;
                        dy2 = _center.y()-tmpCirc._center.y();
                        dy2*=dy2;
                        tmpDist = dx2+dy2;
                        tmpSigma = sqrt( dx2/tmpDist*(_sigmaCenter[0]*_sigmaCenter[0] + tmpCirc._sigmaCenter[0]*tmpCirc._sigmaCenter[0]) +
                                         dy2/tmpDist*(_sigmaCenter[1]*_sigmaCenter[1] + tmpCirc._sigmaCenter[1]*tmpCirc._sigmaCenter[1]) );
                        tmpDist  = sqrt(tmpDist);
                        if (tmpDist<sigmaLevel*tmpSigma) {
                                //cout<<"Compatible by Center"<<endl;
                                return true;
                        }
                }
                return false;
        }

        bool isCompatibleByKrmWith( SimpleCircle2D const& tmpCirc, float sigmaLevel = 1.0 ) {
                float tmpSigma = sqrt( _krmCircFit.covrfd[0] + tmpCirc._krmCircFit.covrfd[0] );
                float tmpDist  = _krmCircFit.rho - tmpCirc._krmCircFit.rho;
                tmpSigma *= sigmaLevel;
                if ( tmpDist>-tmpSigma && tmpDist<tmpSigma) {
                        //cout<<"Compatible by Rho"<<endl;
                        tmpSigma = sqrt( _krmCircFit.covrfd[2] + tmpCirc._krmCircFit.covrfd[2] );
                        tmpDist  = _krmCircFit.phi - tmpCirc._krmCircFit.phi;
                        tmpSigma *= sigmaLevel;
                        if ( tmpDist>-tmpSigma && tmpDist<tmpSigma) {
                                //cout<<"Compatible by Phi"<<endl;
                                tmpSigma = sqrt( _krmCircFit.covrfd[5] + tmpCirc._krmCircFit.covrfd[5] );
                                tmpDist  = _krmCircFit.dca - tmpCirc._krmCircFit.dca;
                                tmpSigma *= sigmaLevel;
                                if ( tmpDist>-tmpSigma && tmpDist<tmpSigma) {
                                        //cout<<"Compatible by Dca"<<endl;
                                        return true;
                                }
                        }
                }
                return false;
        }

        int mergeCirc ( SimpleCircle2D & addCirc, bool &deleteAddCirc, unsigned int minNHitCut = 5) {

                int addedPoints=0;
                bool notPresent;
                std::vector<size_t > hitToRemoveFromAddC;
                size_t iAddCircDHit=0;
                for (std::vector< ptrHitInClosCust >::iterator addCircPoints_it = addCirc._listHitptrs.begin(); addCircPoints_it != addCirc._listHitptrs.end(); ++addCircPoints_it) {
                        notPresent = true;
                        for (std::vector< ptrHitInClosCust >::iterator circPoints_it = _listHitptrs.begin(); circPoints_it != _listHitptrs.end(); ++circPoints_it) {
                                if (addCircPoints_it->getinEventHitID()==circPoints_it->getinEventHitID()) {
                                        notPresent = false;
                                        hitToRemoveFromAddC.push_back(iAddCircDHit);
                                        break;
                                }
                        }
                        if (notPresent) {
                                if (_krmCircFit.mergeHitOfCirc(addCirc._krmCircFit,iAddCircDHit)) {
                                        _listHitptrs.push_back(*addCircPoints_it);
                                        hitToRemoveFromAddC.push_back(iAddCircDHit);
                                        ++addedPoints;
                                }
                        }
                        ++iAddCircDHit;
                }

                float addCircUsefulHits = addCirc._listHitptrs.size() - hitToRemoveFromAddC.size();
                if (addedPoints>0) {
                        float w1, w2, invSumW;
                        float totNhit = _listHitptrs.size() + addCircUsefulHits;
                        w1 = 1.0/(_sigmaRad*_sigmaRad) * (_listHitptrs.size()-addedPoints )/totNhit ;
                        w2 = 1.0/(addCirc._sigmaRad*addCirc._sigmaRad) * addCircUsefulHits/totNhit;
                        invSumW = 1.0/(w1 + w2);
                        _radius = ( (w1*_radius + w2*addCirc._radius)*invSumW );
                        _sigmaRad = sqrt(invSumW);
                        w1 = 1.0/(_sigmaCenter[0]*_sigmaCenter[0]);
                        w2 = 1.0/(addCirc._sigmaCenter[0]*addCirc._sigmaCenter[0]);
                        invSumW = 1.0/(w1 + w2);
                        _center.setX( (w1*_center.x() + w2*addCirc._center.x())*invSumW );
                        _sigmaCenter[0] = sqrt(invSumW);
                        w1 = 1.0/(_sigmaCenter[1]*_sigmaCenter[1]);
                        w2 = 1.0/(addCirc._sigmaCenter[1]*addCirc._sigmaCenter[1]);
                        invSumW = 1.0/(w1 + w2);
                        _center.setY( (w1*_center.y() + w2*addCirc._center.y())*invSumW );
                        _sigmaCenter[1] = sqrt(invSumW);
                }
                if ( addCircUsefulHits > (float)minNHitCut ) {
                        deleteAddCirc = false;
                        for (std::vector<size_t >::iterator hitToRemoveFromAddC_it = hitToRemoveFromAddC.begin(); hitToRemoveFromAddC_it != hitToRemoveFromAddC.end(); ++hitToRemoveFromAddC_it) {
                                addCirc._listHitptrs.erase(addCirc._listHitptrs.begin()+ *hitToRemoveFromAddC_it);
                                //std::vector<circPoint>::iterator removeKrmPnt_it = addCirc._krmCircFit.points.begin()+*hitToRemoveFromAddC_it;
                                //addCirc._krmCircFit.testHit( removeKrmPnt_it->_xx, removeKrmPnt_it->_yy, -removeKrmPnt_it->_errxx, removeKrmPnt_it->_erryy);
                                //addCirc._krmCircFit.points.erase(removeKrmPnt_it);
                                addCirc._krmCircFit.removeHit(*hitToRemoveFromAddC_it);
                        }
                        addCirc._krmCircFit.computeBestCirc();
                } else {
                        deleteAddCirc = true;
                }

                return addedPoints;
        }

        void summCirc ( SimpleCircle2D const & addCirc, int skipFrstNPnts=0 ) {
                if (skipFrstNPnts>addCirc._listHitptrs.size()) {
                        throw cet::exception("RANGE")
                          << "Asked to skip more points than that are present in the added circle"<<endl;

                }
                float w1, w2, invSumW;
                w1 = 1.0/(_sigmaRad*_sigmaRad);
                w2 = 1.0/(addCirc._sigmaRad*addCirc._sigmaRad);
                invSumW = 1.0/(w1 + w2);
                _radius = ( (w1*_radius + w2*addCirc._radius)*invSumW );
                _sigmaRad = sqrt(invSumW);
                w1 = 1.0/(_sigmaCenter[0]*_sigmaCenter[0]);
                w2 = 1.0/(addCirc._sigmaCenter[0]*addCirc._sigmaCenter[0]);
                invSumW = 1.0/(w1 + w2);
                _center.setX( (w1*_center.x() + w2*addCirc._center.x())*invSumW );
                _sigmaCenter[0] = sqrt(invSumW);
                w1 = 1.0/(_sigmaCenter[1]*_sigmaCenter[1]);
                w2 = 1.0/(addCirc._sigmaCenter[1]*addCirc._sigmaCenter[1]);
                invSumW = 1.0/(w1 + w2);
                _center.setY( (w1*_center.y() + w2*addCirc._center.y())*invSumW );
                _sigmaCenter[1] = sqrt(invSumW);

                _listHitptrs.insert(_listHitptrs.end(),(addCirc._listHitptrs.begin()+skipFrstNPnts),addCirc._listHitptrs.end());
                //_krmCircFit.summCirc( addCirc._krmCircFit, skipFrstNPnts );
                _krmCircFit.bestmergeCirc( addCirc._krmCircFit, skipFrstNPnts );
        }

        SimpleCircle2D & operator += ( SimpleCircle2D const & addCirc ) {
                float w1, w2, invSumW;
                w1 = 1.0/(_sigmaRad*_sigmaRad);
                w2 = 1.0/(addCirc._sigmaRad*addCirc._sigmaRad);
                invSumW = 1.0/(w1 + w2);
                _radius = ( (w1*_radius + w2*addCirc._radius)*invSumW );
                _sigmaRad = sqrt(invSumW);
                w1 = 1.0/(_sigmaCenter[0]*_sigmaCenter[0]);
                w2 = 1.0/(addCirc._sigmaCenter[0]*addCirc._sigmaCenter[0]);
                invSumW = 1.0/(w1 + w2);
                _center.setX( (w1*_center.x() + w2*addCirc._center.x())*invSumW );
                _sigmaCenter[0] = sqrt(invSumW);
                w1 = 1.0/(_sigmaCenter[1]*_sigmaCenter[1]);
                w2 = 1.0/(addCirc._sigmaCenter[1]*addCirc._sigmaCenter[1]);
                invSumW = 1.0/(w1 + w2);
                _center.setY( (w1*_center.y() + w2*addCirc._center.y())*invSumW );
                _sigmaCenter[1] = sqrt(invSumW);

                _listHitptrs.insert(_listHitptrs.end(),addCirc._listHitptrs.begin(),addCirc._listHitptrs.end());
                _krmCircFit += addCirc._krmCircFit;

                return *this;
        }

        SimpleCircle2D & operator= ( SimpleCircle2D const& tmpCirc ) {

                _center        = tmpCirc._center;
                _radius        = tmpCirc._radius;
                _sigmaCenter[0]= tmpCirc._sigmaCenter[0];
                _sigmaCenter[1]= tmpCirc._sigmaCenter[1];
                _sigmaRad      = tmpCirc._sigmaRad;
                _listHitptrs   = tmpCirc._listHitptrs;
                _krmCircFit    = tmpCirc._krmCircFit;

        }

        //std::unordered_set< ptrHitInClosCust > _listHitptrs; //list of the pointers to the hits that are associated into the circle
        std::vector< ptrHitInClosCust > _listHitptrs; //list of the pointers to the hits that are associated into the circle

        CLHEP::Hep2Vector _center;
        float _radius;
        float _sigmaCenter[2];
        float _sigmaRad;
        KarimakiCircle _krmCircFit;
};

typedef std::vector<mu2e::SimpleCircle2D> circlesCont;
//typedef std::multimap<size_t, circlesCont::iterator > hitIDCircRel;


struct corssingPoints {
        corssingPoints() :
                _sigmax(0.0),
                _sigmay(0.0),
                _sigmaz(0.0)
        {}

        int const & getRadLayID () const {
                return _hitRef.getRadLayID();
        }

        mu2e::closHitClust const & getHitClust () const {
                return _hitRef.getHitClust();
        }

        int const & getInLayerCellID () const {
                return _hitRef.getinLayerCellID();
        }

        size_t const & getInEventHitID () const {
                return _hitRef.getinEventHitID();
        }

        int const & getCrossRadLayID () const {
                return _crossHitRef.getRadLayID();
        }

        mu2e::closHitClust const & getCrossHitClust () const {
                return _crossHitRef.getHitClust();
        }

        int const & getCrossInLayerCellID () const {
                return _crossHitRef.getinLayerCellID();
        }

        size_t const & getCrossInEventHitID () const {
                return _crossHitRef.getinEventHitID();
        }

        CLHEP::Hep3Vector _pos;
        float _sigmax;
        float _sigmay;
        float _sigmaz;

        ptrHitInClosCust _hitRef;
        ptrHitInClosCust _crossHitRef;
};

typedef std::vector<corssingPoints > points3D;
typedef std::set< std::pair<size_t,size_t>  > hitCrossingList;

typedef std::pair<double,size_t> zHitIdrel;

bool operator <( const zHitIdrel &c1, const zHitIdrel &c2 ) {
        return c1.first < c2.first;
}

typedef std::vector< zHitIdrel > listZHitIdrels;

}  // end namespace mu2e
#endif
