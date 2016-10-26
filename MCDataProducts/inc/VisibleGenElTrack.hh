#ifndef MCDataProducts_VisibleGenElTrack_hh
#define MCDataProducts_VisibleGenElTrack_hh
//
// Data of the Electrons tracks that came from the targets
//
// $Id: VisibleGenElTrack.hh,v 1.4 2013/04/03 22:20:05 tassiell Exp $
// $Author: tassiell $
// $Date: 2013/04/03 22:20:05 $
//
// Original author G. Tassielli
//

// C++ includes.
#include <iostream>
#include <memory>
#include <utility>
#include <map>
#include <vector>
#include <algorithm>

// Mu2e includes.
#include "MCDataProducts/inc/GenElHitData.hh"
//#include "MCDataProducts/inc/SimParticle.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"

// CLHEP includes.
#include "CLHEP/Units/PhysicalConstants.h"

//using namespace art;

namespace mu2e {

typedef art::Ptr<SimParticle> SimParticlePtr;
//typedef std::vector<mu2e::GenElHitData> GenElHitDataCollection;
typedef std::map< art::Ptr<SimParticle>::key_type, mu2e::GenElHitData> GenElHitDataCollection;
typedef std::map<unsigned long, art::Ptr<SimParticle>::key_type, std::less<unsigned long> > GenElHitDCll_OrderByTime;

class VisibleGenElTrack {
public:
        VisibleGenElTrack():_trkKey(0),_trkID(0),_isConversionEl(false),_hitTimeSorted(false),_nTurns(0){
                _hitDataCollection.clear();
                _orederedHitDataColl.clear();
                _turnFrstHit.clear();
                _nHitperTurn.clear();
                _TOFperTurn.clear();
        }
        /*VisibleGenElTrack(SimParticlePtr &trkPtr_ ):
                  _trkPtr(trkPtr_),*/
        VisibleGenElTrack( /*art::ProductID &simPID_,*/ art::RefCore &simRC_, art::Ptr<SimParticle>::key_type &trkKey_, SimParticleCollection::key_type &trkID_, std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> genTrackData_ ):
                //_simPID(simPID_),
                _simRC(simRC_),
                _trkKey(trkKey_),
                _trkID(trkID_),
                _isConversionEl(false),
                _hitTimeSorted(false),
                _nTurns(0)
        {
                //_trkPtr=const_cast<SimParticlePtr&>(trkPtr_);
                //_trkID=_trkPtr->id();
                _genTrackData=genTrackData_;
                _hitDataCollection.clear();
                _orederedHitDataColl.clear();
                _turnFrstHit.clear();
                _nHitperTurn.clear();
                _TOFperTurn.clear();
        }
        ~VisibleGenElTrack(){}

        friend std::ostream& operator<< ( std::ostream& ost,
                        const VisibleGenElTrack& gtd_ ){
                VisibleGenElTrack &gtd = const_cast<VisibleGenElTrack &>(gtd_);
                //ost << "Track id " << gtd._trkPtr.key()<<" vertex "<< gtd._trkPtr->startPosition()<< " Lorentz Vector "<< gtd._trkPtr->startMomentum()<<std::endl;
                ost << "Track id " << gtd._trkID<<" trak key "<<gtd._trkKey<<" vertex "<< gtd._genTrackData.first<< " Lorentz Vector "<< gtd._genTrackData.second<<std::endl;
                ost << "N loops "<<gtd._nTurns<<std::endl;
                for ( unsigned short j=0; j < gtd._nTurns; j++ ) {
                        ost << j+1 <<"-th loop point: "<<gtd._hitDataCollection[ gtd._turnFrstHit[j] ]._hitPoint << " momentum: "<< gtd._hitDataCollection[ gtd._turnFrstHit[j] ]._hitMomentum <<std::endl;
                }
                ost << "Track Hits:"<<std::endl;
                //                  for ( size_t i=0; i < gtd._hitDataCollection.size(); i++ ) {
                //                          ost << gtd._hitDataCollection[i] <<std::endl;
                //                  }
                for ( GenElHitDataCollection::const_iterator hd_it = gtd._hitDataCollection.begin(); hd_it!= gtd._hitDataCollection.end(); ++hd_it ) {
                        ost << hd_it->second <<std::endl;
                }
                return ost;
        }

        void addHitData( StrawHitPtr iHit_, double &mcHitTime_, bool &isOverlapped_, bool &isOvrlpByProton_, bool &isFirst_, CLHEP::Hep3Vector hitPoint_, CLHEP::Hep3Vector hitMomentum_ );/*{
                  GenElHitData hitData_
                  _hitDataCollection.push_back( hitData_ );
          }*/
        void sort();
        bool isNotPresent( size_t iHit_ );
        bool& isConversionEl( );
        void setConversionEl();
        //SimParticlePtr const& getSimPart();
        //art::ProductID& getSimPID();
        art::RefCore& getSimRC();
        SimParticlePtr::key_type& getTrkKey();
        SimParticleCollection::key_type& getTrkID();
        CLHEP::Hep3Vector getTrkVertex();
        CLHEP::HepLorentzVector getTrkLrntzVec();
        size_t getNumOfHit();
        GenElHitData& getHit(/*size_t*/art::Ptr<SimParticle>::key_type ih);
        GenElHitData& getHit(int i_hit);
        GenElHitData& getHitTimeOrder(int i_hit);
        GenElHitData& getFirstHit();
        GenElHitData& getLastHit();
        unsigned short& FindNTurns();
        unsigned short& getNumOfLoops();
        GenElHitData& getithLoopHit(unsigned short ih);
        size_t getithLoopNHit(unsigned short i);
        double getithLoopTOF(unsigned short ih);
        double getInChamberTotPath();
        double getTotPath();

protected:
        //SimParticlePtr _trkPtr;
        //art::ProductID _simPID;
        art::RefCore _simRC;
        SimParticlePtr::key_type _trkKey;
        SimParticleCollection::key_type _trkID;
        std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> _genTrackData;
        bool _isConversionEl;

        GenElHitDataCollection _hitDataCollection;
        GenElHitDCll_OrderByTime _orederedHitDataColl;
        //std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> genTrackData;
        //std::vector hitDataCollection;
        //std::vector::iterator hd_it;
        //typename std::vector::iterator hd_it;
        bool _hitTimeSorted;
        unsigned short _nTurns;
        //std::vector< std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> > turnStartData;
        std::vector< art::Ptr<SimParticle>::key_type/*size_t*/ > _turnFrstHit;
        std::vector< size_t > _nHitperTurn;
        std::vector< double > _TOFperTurn;

        unsigned long rndup(double n);

        GenElHitData _fakeHit;

};

void VisibleGenElTrack::addHitData( StrawHitPtr iHit_, double &mcHitTime_, bool &isOverlapped_, bool &isOvrlpByProton_, bool &isFirst_, CLHEP::Hep3Vector hitPoint_, CLHEP::Hep3Vector hitMomentum_ ) {
        //          _hitDataCollection.push_back(
        //                          GenElHitData( iHit_, mcHitTime_, isOverlapped_, isOvrlpByProton_, isFirst_, hitPoint_, hitMomentum_ )
        //                          );
        if ( _hitDataCollection.insert( GenElHitDataCollection::value_type( iHit_.key(),
                        GenElHitData( iHit_, mcHitTime_, isOverlapped_, isOvrlpByProton_, isFirst_, hitPoint_, hitMomentum_ ) )
                       ).second ) {
                _orederedHitDataColl.insert( GenElHitDCll_OrderByTime::value_type( rndup(mcHitTime_*1000.0), iHit_.key() ) );
        }
}

void VisibleGenElTrack::sort(){
        //std::sort( _hitDataCollection.begin(), _hitDataCollection.end() );
        _hitTimeSorted=true;
}

bool VisibleGenElTrack::isNotPresent( size_t iHit_ ){
        //typename std::vector::iterator hd_it = _hitDataCollection.begin();
        for ( GenElHitDataCollection::iterator hd_it = _hitDataCollection.begin(); hd_it!= _hitDataCollection.end(); ++hd_it ) {
                if(hd_it->first/*_iHit.key()*/==iHit_) return false;
        }
        return  true;
}

bool& VisibleGenElTrack::isConversionEl(){
        return _isConversionEl;
}

void VisibleGenElTrack::setConversionEl(){
        _isConversionEl=true;
}

//SimParticlePtr const& VisibleGenElTrack::getSimPart() {
//       SimParticlePtr const &tmpPtr = _trkPtr;
//       return tmpPtr;
//}

//art::ProductID& VisibleGenElTrack::getSimPID() {
//        return _simPID;
//}

art::RefCore& VisibleGenElTrack::getSimRC() {
        return _simRC;
}

SimParticlePtr::key_type& VisibleGenElTrack::getTrkKey(){
        return _trkKey;
}

SimParticleCollection::key_type& VisibleGenElTrack::getTrkID() {
        //return _trkPtr.key();
        return _trkID;
}

CLHEP::Hep3Vector VisibleGenElTrack::getTrkVertex() {
        //return _trkPtr->startPosition();
        return _genTrackData.first;
}

CLHEP::HepLorentzVector VisibleGenElTrack::getTrkLrntzVec() {
        //return _trkPtr->startMomentum();
        return _genTrackData.second;
}

size_t VisibleGenElTrack::getNumOfHit() {
        return _hitDataCollection.size();
}

GenElHitData& VisibleGenElTrack::getHit(art::Ptr<SimParticle>::key_type/*size_t*/ ih){
        //if ( ih<_hitDataCollection.size() ) return _hitDataCollection[ih];
        GenElHitDataCollection::iterator hd_it=_hitDataCollection.find(ih);
        if ( hd_it!=_hitDataCollection.end() ) return hd_it->second;
        else throw cet::exception("RANGE") << "Element "<<ih<<" out of size of the HitDataCollection "<<_hitDataCollection.size()<<std::endl;
}

GenElHitData& VisibleGenElTrack::getHit(int i_hit){
        if ( ((unsigned int)i_hit)<_hitDataCollection.size() ) {
                GenElHitDataCollection::iterator hd_it=_hitDataCollection.begin();
                for (int i=0; i<i_hit; i++) ++hd_it;
                return hd_it->second;
        }
        else throw cet::exception("RANGE") << "Element "<<i_hit<<" out of size of the HitDataCollection "<<_hitDataCollection.size()<<std::endl;
}

GenElHitData& VisibleGenElTrack::getHitTimeOrder(int i_hit){
        if ( ((unsigned int)i_hit)<_hitDataCollection.size() ) {
                GenElHitDCll_OrderByTime::iterator hd_it=_orederedHitDataColl.begin();
                for (int i=0; i<i_hit; i++) ++hd_it;
                return _hitDataCollection[ hd_it->second ];
        }
        else throw cet::exception("RANGE") << "Element "<<i_hit<<" out of size of the HitDataCollection "<<_hitDataCollection.size()<<std::endl;
}

GenElHitData& VisibleGenElTrack::getFirstHit() {
        return _hitDataCollection[ _orederedHitDataColl.begin()->second ];
}

GenElHitData& VisibleGenElTrack::getLastHit() {
        if ( _hitDataCollection.size()>1) {
                GenElHitDCll_OrderByTime::iterator last_hitto_it = _orederedHitDataColl.end();
                --last_hitto_it;
                return _hitDataCollection[ last_hitto_it->second ];
        } else {
                return _hitDataCollection[ _orederedHitDataColl.begin()->second ];
        }
}

unsigned short& VisibleGenElTrack::FindNTurns(){
        _turnFrstHit.clear();
        _nHitperTurn.clear();
        _TOFperTurn.clear();
        if ( _hitDataCollection.empty() ) {
                std::cerr<<"Not point stored for track to try find loops"<<std::endl;
                _nTurns=0;
                return _nTurns;
        }
        else if ( _hitDataCollection.size() < 10 ) {
                std::cerr<<"Not enough point stored for track to try find loops"<<std::endl;
                _nTurns=1;
		size_t TunrNHit=0;
		_nHitperTurn.push_back(TunrNHit);
                _turnFrstHit.push_back(_orederedHitDataColl.begin()->second);
		_TOFperTurn.push_back(0);
                return _nTurns;
        }
        if ( !_hitTimeSorted ) sort();
        //          for ( size_t ih=1 ; ih<_hitDataCollection.size(); ih++){
        //                  if ( (_hitDataCollection[ih]._mcHitTime - _hitDataCollection[ih-1]._mcHitTime) > 1.0 ) {
        //                          _nTurns++;
        //                          _turnFrstHit.push_back(ih);
        //                  }
        //          }
        _nTurns=1;
        GenElHitDCll_OrderByTime::iterator prev_hitto_it=_orederedHitDataColl.begin();
        _turnFrstHit.push_back(prev_hitto_it->second);
        _TOFperTurn.push_back(_hitDataCollection[ prev_hitto_it->second ]._mcHitTime);

        GenElHitDCll_OrderByTime::iterator hitto_it=prev_hitto_it;
        ++hitto_it;
        size_t TunrNHit=1;
        for ( ; hitto_it!=_orederedHitDataColl.end(); ++hitto_it ){
                //                  if ( (hitto_it->first - prev_hitto_it->first) > 1.0 ) {
                if ( (hitto_it->first - prev_hitto_it->first) > 1200 ) {
                        _TOFperTurn.back()=_hitDataCollection[ prev_hitto_it->second ]._mcHitTime - _TOFperTurn.back();
                        _nHitperTurn.push_back(TunrNHit);
                        _nTurns++;
                        _turnFrstHit.push_back(hitto_it->second);
                        _TOFperTurn.push_back(_hitDataCollection[ hitto_it->second ]._mcHitTime);
                        TunrNHit=1;
                }
                prev_hitto_it=hitto_it;
                TunrNHit++;
        }
        --hitto_it;
        _TOFperTurn.back()=_hitDataCollection[ hitto_it->second ]._mcHitTime-_TOFperTurn.back();
        _nHitperTurn.push_back(TunrNHit);
        return _nTurns;
}

unsigned short& VisibleGenElTrack::getNumOfLoops(){
        if (_nTurns==0) { FindNTurns(); }
        return _nTurns;
}

GenElHitData& VisibleGenElTrack::getithLoopHit(unsigned short ih){
        //if (ih>=_nTurns) throw cet::exception("RANGE") << "Requested Loop "<<ih<<" -th not found for the track"<<std::endl;
        if (ih>=_nTurns) {
                std::cerr << "Requested Loop "<<ih<<" -th not found for the track"<<std::endl;
                return _fakeHit;
        }
        return _hitDataCollection[ _turnFrstHit[ih] ];
}

size_t VisibleGenElTrack::getithLoopNHit(unsigned short ih){
        if (ih>=_nTurns) return 0;
        return _nHitperTurn.at(ih);
}

double VisibleGenElTrack::getithLoopTOF(unsigned short ih){
        if (ih>=_nTurns) return 0.0;
        return _TOFperTurn.at(ih);
}

double VisibleGenElTrack::getInChamberTotPath(){
        if (_nTurns==0) { FindNTurns(); }
        double TotalPath(0.0);
        for (std::vector< double >::iterator TOFperTurn_it = _TOFperTurn.begin(); TOFperTurn_it != _TOFperTurn.end(); ++TOFperTurn_it) {
                TotalPath+=*TOFperTurn_it;
        }
        TotalPath*=CLHEP::c_light;
        return TotalPath;
}

double VisibleGenElTrack::getTotPath(){
        if (_nTurns==0) { FindNTurns(); }
        double TotalPath(0.0);
        if ( _hitDataCollection.size()>1) {
                GenElHitDCll_OrderByTime::iterator last_hitto_it = _orederedHitDataColl.end();
                --last_hitto_it;
                TotalPath = _hitDataCollection[ last_hitto_it->second ]._mcHitTime - _hitDataCollection[ _orederedHitDataColl.begin()->second ]._mcHitTime;
                TotalPath*=CLHEP::c_light;
        }
        return TotalPath;
}

unsigned long VisibleGenElTrack::rndup(double n)//round up a float type and show one decimal place
{
        double t;
        t=n-floor(n);
        if (t>=0.5)
        {
                n*=10.00000;//where n is the multi-decimal float
                n=ceil(n);
                n/=10.00000;
        }
        else
        {
                n*=10.00000;//where n is the multi-decimal float
                n=floor(n);
                n/=10.00000;
        }
        return (unsigned long)n;
}


}

#endif /* MCDataProducts_VisibleGenElTrack_hh */
