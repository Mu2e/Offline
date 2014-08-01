//
// General utilities for the calorimeter's studies
//
// $Id: CaloClusterUtilities.hh,v 1.10 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Original author G. Pezzullo & G. Tassielli & G. Onorato
//

//art includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

//Mu2e includes
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CaloCluster/inc/CaloClusterTools.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "CLHEP/Vector/ThreeVector.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/Exception.h"
#include "art/Framework/Principal/Event.h"

// Mu2e includes

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "BaBar/BbrGeom/include/BbrVectorErr.hh"
#include "BaBar/BbrGeom/include/BbrPointErr.hh"

//c++ includes
#include <string>

//Root includes
#include "TMath.h"

using namespace std;

namespace mu2e {


class MCCaloUtilities {

public:

        MCCaloUtilities();

        ~MCCaloUtilities();

        void setTrackAndRO(const art::Event & event,
                        std::string const &_g4ModuleLabel,
                        SimParticleCollection::key_type track,
                        unsigned RO);

        void printOutCaloInfo();

        bool fromOutside();

        bool primary();

        bool generated();

        int startingVane();

        int getStartingVane(CLHEP::Hep3Vector origin);

        int localVane();

private:

        unsigned _localRO;
        unsigned _localCrystal;
        unsigned _localVane;
        int _startingVane;
        bool _fromOutside, _primary, _generated;

};

std::string & TOUpper(std::string &in);

double cry(double val);

double indexToCoor(double ind);
double indexToCoor(int ind);

void w_correction_0( double& clCOGw,  double& clCOGwErr,  int& clCryEnergyMaxColumn);

void v_correction_0( float& extrapolThetaV,  double& clCOGv,  double& clCOGwErr);

void v_correction_0( double& extrapolThetaV,  double& clCOGv,  double& clCOGwErr);

void w_correction_1(double& clCOGw,double& clCOGwErr, int& wSize);

void v_correction_1(int& clCryEnergyMaxRow,  double& clCOGv,  double& clCOGvErr);

void cog_correction_0(CaloCluster &cluster);


//define a map which key is the index of the row "R" of the vane, and the contained object is an other map which key is the column index "Z" and also contain a vector.
//The vector contains pairs of "CaloCrystalHit" and a index which is the position of the "CaloCrystalHit" in the vector "CaloCrystalHitCollection" generated in the event
typedef std::map<unsigned int, std::map<unsigned int, std::vector<std::pair<CaloCrystalHit *, size_t > > >  > MatrixCaloHit;

//define a map which key is the vane's index and contains object of type "CaloCrystalHit". In that way we have a complete topology description of the calorimeter
typedef std::map<unsigned int, MatrixCaloHit> VanesMap;


//define the object in which we store a single cluster. the key is the row index, and the pair contains the column index and the position of the CaloCrystalHit in the vector
//stored in the container "vanesMap"
typedef std::multimap<unsigned int, std::pair<unsigned int, unsigned int> > ClusterData;//row, cloumn, hitId

struct DistVector{
        std::vector<double> _vec;

        void push_back(double id){
                _vec.push_back(id);
        }
        size_t size(){
                return _vec.size();
        }
        bool find( double t){
                bool res = false;

                unsigned int size = _vec.size();
                if(size!=0){
                        unsigned int cont = 0;
                        while( cont!=size){
                                if(_vec[cont] == t) {
                                        res = true;
                                }
                                ++cont;
                        }
                }
                return res;
        }

        size_t indexMinCont(){
                size_t pos = 0;
                if(_vec.size() == 0) return pos;

                double tmpMin = _vec.at(0);

                for(size_t i =0; i< _vec.size(); ++i){
                        if(_vec[i] < tmpMin) {
                                tmpMin = _vec[i];
                                pos = i;
                        }
                }

                return pos;
        }

        //Accessors
        double &operator[](size_t n){
                return (this->_vec.at(n));
        }
        double &at(size_t &n){
                return (this->_vec.at(n));
        }


        DistVector & operator=(const DistVector other){
                for(size_t i =0; i<other._vec.size(); ++i){
                        _vec.push_back( other._vec.at(i) );
                }
                return *this;
        }

        DistVector(){};

};

struct trkIdVector{
        std::vector<unsigned int> _vec;

        void push_back(int id){
                _vec.push_back((unsigned int)id);
        }

        size_t size(){
                return _vec.size();
        }
        bool find( unsigned int t){
                bool res = false;
                if(_vec.size() == 0) return res;

                for(size_t cont=0; cont<_vec.size(); ++cont){
                        if(_vec[cont] == t) {
                                res = true;
                        }
                }

                return res;
        }

        void remove(unsigned int index){
                if( (_vec.size() == 0) || !(this->find(index)) ) return ;

                trkIdVector resVec;
                for(size_t j=0; j< _vec.size(); ++j){
                        if(!(this->find(index)) ){
                                resVec.push_back(_vec.at(j));
                        }
                }
                this->_vec.clear();
                *this = resVec;
        }

        //Accessors
        unsigned int &operator[](size_t n){
                return (this->_vec.at(n));
        }
        unsigned int &at(size_t &n){
                return (this->_vec.at(n));
        }


        trkIdVector & operator=(const trkIdVector other){
                for(size_t i =0; i<other._vec.size(); ++i){
                        _vec.push_back( other._vec.at(i) );
                }
                return *this;
        }
        trkIdVector(){}

        void print(std::ostream& os) const{
                os << "trkIdVector print out :"<<std::endl;
                for(size_t i=0; i<_vec.size(); ++i){
                        os << "vec("<<i<<") = "<<_vec.at(i)<<std::endl;
                }
        }

};

struct elecData{
  int _vaneId, _isConv, _pdgId, _isGen;
        double _impTime, _impTimeErr;
        double _genTime, _genTimeErr;
        int _index;
        CLHEP::Hep3Vector _impPos;
        BbrPointErr _impPosErr;
        CLHEP::Hep3Vector _impMom3Vec;
        BbrVectorErr _impMom3VecErr;
        double  _pathLenght
        , _t0
        , _pathLenghtErr;
        double _t0Err;
        CLHEP::Hep3Vector _t0Momentum;
        BbrVectorErr _t0MomentumErr;
        CLHEP::HepLorentzVector _genMom;
        double _fitConsistency;

        bool operator<( const elecData other) const{
                return ( _impTime< other._impTime);
        }

        bool operator>( const elecData other) const{
                return ( _impTime> other._impTime);
        }

        bool operator==( const elecData other) const{
                if( _vaneId == other._vaneId &&
		    _isConv == other._isConv &&
		    _pdgId == other._pdgId &&
		    _isGen == other._isGen &&
                                _index   == other._index &&
                                _impTime == other._impTime &&
                                _impTimeErr == other._impTimeErr &&
		    _genTime == other._genTime &&
                                _genTimeErr == other._genTimeErr &&
                                _impPos == other._impPos &&
                                _impPosErr == other._impPosErr &&
                                _impMom3Vec == other._impMom3Vec &&
                                _impMom3VecErr == other._impMom3VecErr &&
                                _t0Momentum == other._t0Momentum &&
                                _t0MomentumErr == other._t0MomentumErr &&
                                _pathLenght == other._pathLenght &&
                                _pathLenghtErr == other._pathLenghtErr &&
                                _t0Err == other._t0Err &&
                                _t0 == other._t0 &&
                                _genMom == other._genMom &&
                                _fitConsistency== other._fitConsistency){
                        return true;
                }else {
                        return false;
                }
        }
        elecData & operator=(const elecData& other) {
                _vaneId         = other._vaneId;
		_isConv         = other._isConv;
		_pdgId         = other._pdgId;
		_isGen         = other._isGen;
                _index          = other._index;
                _impTime        = other._impTime;
                _impTimeErr     = other._impTimeErr;
		_genTime        = other._genTime;
                _genTimeErr     = other._genTimeErr;
		_impPos         = other._impPos;
                _impPosErr      = other._impPosErr;
                _impMom3Vec     = other._impMom3Vec;
                _impMom3VecErr  = other._impMom3VecErr;
                _t0Momentum       = other._t0Momentum;
                _t0MomentumErr      = other._t0MomentumErr;
                _pathLenght     = other._pathLenght;
                _pathLenghtErr  = other._pathLenghtErr;
                _t0             = other._t0;
                _t0Err             = other._t0Err;
                _genMom             = other._genMom;
                _fitConsistency    = other._fitConsistency;

                return *this;
        }
        elecData():
	  _isConv(0),
	  _pdgId(0),
	  _isGen(0),
                _impTime(1e10),
                _impTimeErr(0.0),
	  _genTime(0.0),
	  _genTimeErr(0.0),
                _index(0),
                _pathLenght(0.0),
                _t0(0.0),
                _pathLenghtErr(0.0),
                _t0Err(0.0),
                _fitConsistency(0.0){
        }

        void print(std::ostream& os) const{
                os << "vaneId = "<<_vaneId <<
		  ", isConv = "<< _isConv <<
		  ", pdgId = "<< _pdgId <<
		  ", isGen = "<< _isGen <<
		  
                                ", index = "<<_index <<
                                ", impTime = "<<_impTime << " [ns]"<<
                                ", impTimeErr = "<<_impTimeErr << " [ns]"<<
		  ",genTime = "<<_genTime << " [ns]"<<
                                ", genTimeErr = "<<_genTimeErr << " [ns]"<<
                                ", impPos = "<<_impPos << " [mm]"<<
                                ", pathLength = "<<_pathLenght << " [mm]"<<
                                ", pathLengthErr = "<<_pathLenghtErr << " [mm]"<<
                                ", t0 = "<<_t0 << " [ns]"<<
                                ", t0Err = "<<_t0Err << " [ns]"<<
                                ", t0Momentum = "<<_t0Momentum << " [ns]"<<
                                ", t0MomentumErr = "<<_t0MomentumErr << " [ns]"<<
                                ", impPosErr = "<<_impPosErr << " [mm]"<< std::endl;
                os << "impMom3Vec = "<<_impMom3Vec << " [MeV]"<<
                                ", impMom3VecErr = "<<_impMom3VecErr << " [MeV]"<<
                                ", energy = "<<_impMom3Vec.mag() << " [MeV]"<<
                                ", generation energy = "<<_genMom.e() << " [MeV]"<<
                                ", fitConsistency = "<<_fitConsistency <<std::endl;

        }
};


struct elecDataVector{
        std::vector< elecData > _vec;

        void push_back(elecData &id){
                _vec.push_back(id);
        }
        size_t size(){
                return _vec.size();
        }

        //returns the vector which contains the list of the vaneId intersected
        trkIdVector vaneList()const{
                trkIdVector res;

                if(_vec.size() == 0)return res;

                for(size_t j=0; j<_vec.size(); ++j){
                        unsigned int vaneIndex = (unsigned int) _vec.at(j)._vaneId;
                        if(!res.find( vaneIndex) ) {
                                res.push_back(vaneIndex);
                        }
                }
                return res;
        }

        //returns a std::vector<> which contains the indexes which belong to the same vaneId
        trkIdVector findVane( int vane){
                trkIdVector res;

                if(_vec.size() == 0)return res;

                for(size_t j=0; j<_vec.size(); ++j){
                        if(_vec.at(j)._vaneId == vane) {
                                res.push_back(j);
                        }
                }
                return res;
        }

        bool findElecData(elecData &id, size_t &index){
                bool res = false;
                if(_vec.size()==0) return res;
                size_t count = 0;
                while( ( count != _vec.size() ) && !res){
                        if(_vec.at(count) == id){
                                res = true;
                        }else{
                                ++count;
                        }
                }
                index = count;
                return res;
        }

        void remove(elecData &id){
                bool res = false;
                size_t count = 0;
                while( ( count != _vec.size() ) && !res){
                        if(_vec.at(count) == id){
                                res = true;
                        }else{
                                ++count;
                        }
                }
                size_t index = count;
                if( (_vec.size() == 0) ) return ;

                elecDataVector resVec;
                for(size_t j=0; j< _vec.size(); ++j){
                        if( j != index ){
                                resVec.push_back(_vec.at(j));
                        }
                }
                this->_vec.clear();
                *this = resVec;
        }

        void remove(size_t &index){

                if( (_vec.size() == 0) ) return ;

                elecDataVector resVec;
                for(size_t j=0; j< _vec.size(); ++j){
                        if( j != index ){
                                resVec.push_back(_vec.at(j));
                        }
                }
                this->_vec.clear();
                *this = resVec;
        }
        //Accessors
        elecData &operator[](size_t n){
                return (this->_vec.at(n));
        }
        elecData &at(size_t &n){
                return (this->_vec.at(n));
        }
        int                         vaneId(size_t &n)const{return _vec.at(n)._vaneId;}
        int                         isConv(size_t &n)const{return _vec.at(n)._isConv;}
        int                         pdgId(size_t &n)const{return _vec.at(n)._pdgId;}
        int                         isGen(size_t &n)const{return _vec.at(n)._isGen;}
        int                        index(size_t &n)const{return _vec.at(n)._index;}
        double                     impTime(size_t &n)const{return _vec.at(n)._impTime;}
        double                        impTimeErr(size_t &n)const{return _vec.at(n)._impTimeErr;}
        CLHEP::Hep3Vector                impPos(size_t& n)const{return _vec.at(n)._impPos;}
        BbrPointErr            impPosErr(size_t& n)const{return _vec.at(n)._impPosErr;}
        CLHEP::Hep3Vector        impMom3Vec(size_t& n)const{return _vec.at(n)._impMom3Vec;}
        BbrVectorErr           impMom3VecErr(size_t& n)const{return _vec.at(n)._impMom3VecErr;}
        double                        t0(size_t &n)const{return _vec.at(n)._t0;}
        double                        t0Err(size_t &n)const{return _vec.at(n)._t0Err;}
        CLHEP::Hep3Vector                   t0Momentum(size_t &n)const{return _vec.at(n)._t0Momentum;}
        BbrVectorErr               t0MomentumErr(size_t &n)const{return _vec.at(n)._t0MomentumErr;}
        double                   pathLenght(size_t &n)const{return _vec.at(n)._pathLenght;}
        double                   pathLenghtErr(size_t &n)const{return _vec.at(n)._pathLenghtErr;}


        elecData &operator[](unsigned int n){
                return (this->_vec.at(n));
        }
        elecData &at(unsigned int &n){
                return (this->_vec.at(n));
        }
  int                          vaneId(unsigned int &n)const{return _vec.at(n)._vaneId;}
  int                          isConv(unsigned int &n)const{return _vec.at(n)._isConv;}
  int                         isGen(unsigned int &n)const{return _vec.at(n)._isGen;}
        int                         pdgId(unsigned int &n)const{return _vec.at(n)._pdgId;}
        int                            index(unsigned int &n)const{return _vec.at(n)._index;}
        double                        impTime(unsigned int &n)const{return _vec.at(n)._impTime;}
        double                        impTimeErr(unsigned int &n)const{return _vec.at(n)._impTimeErr;}
        CLHEP::Hep3Vector                impPos(unsigned int& n)const{return _vec.at(n)._impPos;}
        BbrPointErr                    impPosErr(unsigned int& n)const{return _vec.at(n)._impPosErr;}
        CLHEP::Hep3Vector            impMom3Vec(unsigned int& n)const{return _vec.at(n)._impMom3Vec;}
        BbrVectorErr                  impMom3VecErr(unsigned int& n)const{return _vec.at(n)._impMom3VecErr;}
        double                          t0(unsigned int &n)const{return _vec.at(n)._t0;}
        double                        t0Err(unsigned int &n)const{return _vec.at(n)._t0Err;}
        CLHEP::Hep3Vector                    t0Momentum(unsigned int &n)const{return _vec.at(n)._t0Momentum;}
        BbrVectorErr               t0MomentumErr(unsigned int &n)const{return _vec.at(n)._t0MomentumErr;}
        double                  pathLenght(unsigned int &n)const{return _vec.at(n)._pathLenght;}
        double                   pathLenghtErr(unsigned int &n)const{return _vec.at(n)._pathLenghtErr;}


        elecDataVector & operator=(const elecDataVector other){
                for(size_t i =0; i<other._vec.size(); ++i){
                        _vec.push_back( other._vec.at(i) );
                }
                return *this;
        }

        elecDataVector(){};

};


//the key is track Id and the second key is the track generation Id
typedef std::map<unsigned int,elecDataVector > ElecMap;

//---------------------------------------------------------------

struct clusterData{
        int _trkId;
        double _cluTime;
        CLHEP::Hep3Vector _cluPos, _cluPosErr;
        double _cluEnergy, _cluEnergyErr;

        bool operator<( const clusterData other) const{
                return ( _cluTime< other._cluTime);
        }

        bool operator>( const clusterData other) const{
                return ( _cluTime> other._cluTime);
        }

        bool operator==( const clusterData other) const{
                if( _trkId == other._trkId &&
                                _cluTime == other._cluTime &&
                                _cluPos == other._cluPos &&
                                _cluPosErr == other._cluPosErr &&
                                _cluEnergy == other._cluEnergy &&
                                _cluEnergyErr == other._cluEnergyErr ){
                        return true;
                }else {
                        return false;
                }
        }
        clusterData & operator=(const clusterData& other) {
                _trkId         = other._trkId;
                _cluTime        = other._cluTime;
                _cluPos         = other._cluPos;
                _cluPosErr      = other._cluPosErr;
                _cluEnergy     = other._cluEnergy;
                _cluEnergyErr  = other._cluEnergyErr;

                return *this;
        }
        clusterData():
                _cluTime(1e10){
        }

        void print(std::ostream& os) const{
                os << "vaneId = "<<_trkId <<
                                ", cluTime = "<<_cluTime << " [ns]"<<
                                ", cluPos = "<<_cluPos << " [mm]"<<
                                ", cluPosErr = "<<_cluPosErr << " [mm]"<< std::endl;
                os << "cluMom3Vec = "<<_cluEnergy << " [MeV]"<<
                                ", cluMom3VecErr = "<<_cluEnergyErr << " [MeV]"<<
                                ", energy = "<<_cluEnergy << " [MeV]"<< std::endl;
        }
};


struct clusterDataVector{
        std::vector< clusterData > _vec;

        void push_back(clusterData &id){
                _vec.push_back(id);
        }
        size_t size(){
                return _vec.size();
        }

        //returns the vector which contains the list of the trkId intersected
        trkIdVector trkIdList()const{
                trkIdVector res;

                if(_vec.size() == 0)return res;

                for(size_t j=0; j<_vec.size(); ++j){
                        unsigned int trkIdIndex = (unsigned int) _vec.at(j)._trkId;
                        if(!res.find( trkIdIndex) ) {
                                res.push_back(trkIdIndex);
                        }
                }
                return res;
        }

        //returns a std::vector<> which contains the indexes which belong to the same trkId
        trkIdVector findtrkId( int trkId){
                trkIdVector res;

                if(_vec.size() == 0)return res;

                for(size_t j=0; j<_vec.size(); ++j){
                        if(_vec.at(j)._trkId == trkId) {
                                res.push_back(j);
                        }
                }
                return res;
        }

        bool findclusterData(clusterData &id, size_t &index){
                bool res = false;
                if(_vec.size()==0) return res;
                size_t count = 0;
                while( ( count != _vec.size() ) && !res){
                        if(_vec.at(count) == id){
                                res = true;
                        }else{
                                ++count;
                        }
                }
                index = count;
                return res;
        }

        void remove(clusterData &id){
                bool res = false;
                size_t count = 0;
                while( ( count != _vec.size() ) && !res){
                        if(_vec.at(count) == id){
                                res = true;
                        }else{
                                ++count;
                        }
                }
                size_t index = count;
                if( (_vec.size() == 0) ) return ;

                clusterDataVector resVec;
                for(size_t j=0; j< _vec.size(); ++j){
                        if( j != index ){
                                resVec.push_back(_vec.at(j));
                        }
                }
                this->_vec.clear();
                *this = resVec;
        }

        void remove(size_t &index){

                if( (_vec.size() == 0) ) return ;

                clusterDataVector resVec;
                for(size_t j=0; j< _vec.size(); ++j){
                        if( j != index ){
                                resVec.push_back(_vec.at(j));
                        }
                }
                this->_vec.clear();
                *this = resVec;
        }
        //Accessors
        clusterData &operator[](size_t n){
                return (this->_vec.at(n));
        }
        clusterData &at(size_t &n){
                return (this->_vec.at(n));
        }
        int trkId(size_t &n)const{return _vec.at(n)._trkId;}
        double cluTime(size_t &n)const{return _vec.at(n)._cluTime;}
        CLHEP::Hep3Vector cluPos(size_t& n)const{return _vec.at(n)._cluPos;}
        CLHEP::Hep3Vector cluPosErr(size_t& n)const{return _vec.at(n)._cluPosErr;}
        double cluEnergy(size_t& n)const{return _vec.at(n)._cluEnergy;}
        double cluEnergyErr(size_t& n)const{return _vec.at(n)._cluEnergyErr;}

        clusterData &operator[](unsigned int n){
                return (this->_vec.at(n));
        }
        clusterData &at(unsigned int &n){
                return (this->_vec.at(n));
        }
        int trkId(unsigned int &n)const{return _vec.at(n)._trkId;}
        double cluTime(unsigned int &n)const{return _vec.at(n)._cluTime;}
        CLHEP::Hep3Vector cluPos(unsigned int& n)const{return _vec.at(n)._cluPos;}
        CLHEP::Hep3Vector cluPosErr(unsigned int& n)const{return _vec.at(n)._cluPosErr;}
        double cluEnergy(unsigned int& n)const{return _vec.at(n)._cluEnergy;}
        double cluEnergyErr(unsigned int& n)const{return _vec.at(n)._cluEnergyErr;}

        clusterDataVector & operator=(const clusterDataVector other){
                for(size_t i =0; i<other._vec.size(); ++i){
                        _vec.push_back( other._vec.at(i) );
                }
                return *this;
        }

        clusterDataVector(){};

};


//the key is track Id and the second key is the vane index
typedef std::map<unsigned int,clusterDataVector > ClusterDataMap;

struct ClusterMap{
        int _cluSize;
        int _COGcrySize;
        std::vector<double> _rowVec;
        std::vector<double> _columnVec;
        std::vector<double> _cryEdepVec;
        std::vector<double> _COGrowVec;
        std::vector<double> _COGcolumnVec;
        std::vector<double> _timeVec;
        int _vaneId;
        int _cluCogRow;
        int _cluCogColumn;
        int _cryEnergydepMaxRow;
        int _cryEnergydepMaxColumn;
        CLHEP::Hep3Vector _cluCOG;
        float   _time;
        float   _showerDir;
        float   _errShowerDir;

        void print(std::ostream& os) const{
                os << "vaneId = "<<_vaneId <<
                                ", cluSize = "<<_cluSize << " [# crystals]"<<
                                ", time = "<<_time << " [ns]"<<
                                ", cluEnergy = "<< this->clusterEnergy()<<" [MeV]"<<
                                ", COGcrySize = "<<_COGcrySize << " [# crystals]"<<endl<<
                                ", cluCogRow = "<<_cluCogRow <<
                                ", cluCogColumn = "<<_cluCogColumn <<
                                ", cryEnergydepMaxRow = "<<_cryEnergydepMaxRow <<
                                ", cryEnergydepMaxColumn = "<<_cryEnergydepMaxColumn <<endl<<
                                ", showerDir = "<<_showerDir <<
                                ", errShowerDir = "<<_errShowerDir <<std::endl;

                for(size_t i=0; i<_rowVec.size(); ++i){
                        os<< "rowVec["<<i<<"] = "<<_rowVec[i]<<std::endl;
                }
                for(size_t i=0; i<_rowVec.size(); ++i){
                        os<< "columnVec["<<i<<"] = "<<_columnVec[i]<<std::endl;
                }
                for(size_t i=0; i<_rowVec.size(); ++i){
                        os<< "cryEdepVec["<<i<<"] = "<<_cryEdepVec[i]<<std::endl;
                }
                for(size_t i=0; i<_rowVec.size(); ++i){
                        os<< "COGrowVec["<<i<<"] = "<<_COGrowVec[i]<<std::endl;
                }
                for(size_t i=0; i<_rowVec.size(); ++i){
                        os<< "COGcolumnVec["<<i<<"] = "<<_COGcolumnVec[i]<<std::endl;
                }
		for(size_t i=0; i<_rowVec.size(); ++i){
                        os<< "timeVec["<<i<<"] = "<<_timeVec[i]<<std::endl;
                }
        }

        bool operator==( const ClusterMap other) const{
                if(_cluSize             == other._cluSize &&
                                _COGcrySize             == other._COGcrySize &&
                                _time                   == other._time &&
                                _vaneId                 == other._vaneId &&
                                _cluCogRow              == other._cluCogRow &&
                                _cluCogColumn           == other._cluCogColumn &&
                                _cluCOG                 == other._cluCOG &&
                                _rowVec                 == other._rowVec &&
                                _columnVec              == other._columnVec &&
                                _cryEdepVec             == other._cryEdepVec &&
                                _COGrowVec              == other._COGcolumnVec &&
                                _COGcolumnVec           == other._COGcolumnVec &&
		                _timeVec                == other._timeVec &&
                                _showerDir              == other._showerDir &&
                                _errShowerDir           == other._errShowerDir &&
                                _cryEnergydepMaxRow     == other._cryEnergydepMaxRow &&
                                _cryEnergydepMaxColumn  == other._cryEnergydepMaxColumn){
                        return true;
                }else {
                        return false;
                }
        }
        ClusterMap & operator=(const ClusterMap& other) {
                _cluSize                = other._cluSize;
                _COGcrySize             = other._COGcrySize;
                _time                   = other._time ;
                _vaneId                 = other._vaneId;
                _cluCogRow              = other._cluCogRow;
                _cluCogColumn           = other._cluCogColumn;
                _cluCOG                 = other._cluCOG;
                _rowVec                 = other._rowVec;
                _columnVec              = other._columnVec;
                _cryEdepVec             = other._cryEdepVec;
                _COGrowVec              = other._COGcolumnVec;
                _COGcolumnVec           = other._COGcolumnVec;
		_timeVec                = other._timeVec;
                _showerDir              = other._showerDir;
                _errShowerDir           = other._errShowerDir;
                _cryEnergydepMaxRow     = other._cryEnergydepMaxRow;
                _cryEnergydepMaxColumn  = other._cryEnergydepMaxColumn;



                return *this;
        }
        ClusterMap():
                _cluSize(0),
                _cluCogRow(0.0),
                _cluCogColumn(0.0),
                _cryEnergydepMaxRow(0),
                _cryEnergydepMaxColumn(0),
                _time(0.0){
        }
        ClusterMap(CaloCluster &cluster):
                _cluSize(cluster.size()),
                _COGcrySize(cluster.size()),
                _vaneId(cluster.vaneId()),
                _cluCogRow(cluster.cogRow()),
                _cluCogColumn(cluster.cogColumn()),
                _cryEnergydepMaxRow(CaloClusterTools(cluster).cryEnergydepMaxRow() ),
                _cryEnergydepMaxColumn(CaloClusterTools(cluster).cryEnergydepMaxColumn() ),
                _cluCOG(cluster.cog3Vector()),
                _time(cluster.time()),
                _showerDir(CaloClusterTools(cluster).showerDir()),
                _errShowerDir(CaloClusterTools(cluster).errShowerDir()){
                GeomHandle<VaneCalorimeter> cg;
                for( CaloCrystalHitPtrVector::const_iterator itCD = cluster.caloCrystalHitsPtrVector().begin(); itCD != cluster.caloCrystalHitsPtrVector().end(); ++itCD){
		  std::vector<art::Ptr<CaloHit> > const& ROIds = (*itCD)->readouts();
		  if (ROIds.size() <= 0) {
//-----------------------------------------------------------------------------
// 2013-04-15 P.Murat: work around some error, diagnosing which seems difficult at this point
//-----------------------------------------------------------------------------
		    printf(">>> ERROR in CaloClusterUtilities::ClusterMap(CaloCluster &cluster): empty ROIds, initialization incomplete\n");
		  }
		  else {

		    CaloHit const& thehit = *ROIds.at(0);
		    int roid = thehit.id();
		    int tmpRow = cg->crystalRByRO( roid) ;
		    int tmpColumn = cg->crystalZByRO(roid) ;
		    _rowVec.push_back(tmpRow);
		    _columnVec.push_back(tmpColumn);
		    _cryEdepVec.push_back((*itCD)->energyDep());
		    _COGrowVec.push_back(tmpRow);
		    _COGcolumnVec.push_back(tmpColumn);
		    _timeVec.push_back((*itCD)->time());
		  }
                }
        }

        //Setting
        void   setUcog(double newU){_cluCOG.setX(newU);}
        void   setVcog(double newV){
                GeomHandle<VaneCalorimeter> cg;
                _cluCOG.setY(newV);
                int tmpRow=0;
                tmpRow = (int)(newV/cg->caloGeomInfo().crystalHalfTrans());
                _cluCogRow = tmpRow;
        }
        void   setWcog(double newW){
                GeomHandle<VaneCalorimeter> cg;
                _cluCOG.setZ(newW);
                int tmpColumn=0;
                tmpColumn = (int)(newW/cg->caloGeomInfo().crystalHalfTrans());
                _cluCogColumn = tmpColumn;
        }

        //Accessors
        double Ucog() const{return _cluCOG.x();}
        double Vcog() const{return _cluCOG.y();}
        double Wcog() const{return _cluCOG.z();}
        double clusterEnergy() const{
                double E=0.0;
                for(size_t i=0; i<_cryEdepVec.size(); ++i){
                        E += _cryEdepVec.at(i);
                }
                return E;
        }
};


struct ClusterMapVector{
        std::vector<ClusterMap> _vec;

        void push_back(ClusterMap &id){
                _vec.push_back(id);
        }
        size_t size(){
                return _vec.size();
        }

        //returns the vector which contains the list of the vaneId intersected
        trkIdVector vaneList()const{
                trkIdVector res;

                if(_vec.size() == 0)return res;

                for(size_t j=0; j<_vec.size(); ++j){
                        unsigned int vaneIndex = (unsigned int) _vec.at(j)._vaneId;
                        if(!res.find( vaneIndex) ) {
                                res.push_back(vaneIndex);
                        }
                }
                return res;
        }

        //returns a std::vector<> which contains the indexes which belong to the same vaneId
        trkIdVector findVane( int vane){
                trkIdVector res;

                if(_vec.size() == 0)return res;

                for(size_t j=0; j<_vec.size(); ++j){
                        if(_vec.at(j)._vaneId == vane) {
                                res.push_back(j);
                        }
                }
                return res;
        }

        bool searchVane(int vane){
                bool res = true;

                size_t count = 0;
                while(!res && count!=_vec.size()){
                        if(_vec.at(count)._vaneId == vane){
                                res = true;
                        }else{
                                ++count;
                        }
                }
                return res;
        }

        bool findClusterMap(ClusterMap &id, size_t &index){
                bool res = false;
                if(_vec.size()==0) return res;
                size_t count = 0;
                while( ( count != _vec.size() ) && !res){
                        if(_vec.at(count) == id){
                                res = true;
                        }else{
                                ++count;
                        }
                }
                index = count;
                return res;
        }

        void remove(ClusterMap &id){
                bool res = false;
                size_t count = 0;
                while( ( count != _vec.size() ) && !res){
                        if(_vec.at(count) == id){
                                res = true;
                        }else{
                                ++count;
                        }
                }
                size_t index = count;
                if( (_vec.size() == 0) ) return ;

                ClusterMapVector resVec;
                for(size_t j=0; j< _vec.size(); ++j){
                        if( j != index ){
                                resVec.push_back(_vec.at(j));
                        }
                }
                this->_vec.clear();
                *this = resVec;
        }

        void remove(size_t &index){

                if( (_vec.size() == 0) ) return ;

                ClusterMapVector resVec;
                for(size_t j=0; j< _vec.size(); ++j){
                        if( j != index ){
                                resVec.push_back(_vec.at(j));
                        }
                }
                this->_vec.clear();
                *this = resVec;
        }

        //Accessors
        ClusterMap &operator[](size_t n){
                return (this->_vec.at(n));
        }
        ClusterMap &at(size_t &n){
                return (this->_vec.at(n));
        }

        int                                     clusterSize(size_t &n){return _vec.at(n)._cluSize;}
        int                                     COGcrystalSize(size_t &n){return _vec.at(n)._COGcrySize;}
        std::vector<double>                     rowVector(size_t &n){return _vec.at(n)._rowVec;}
        std::vector<double>                     columnVector(size_t &n){return _vec.at(n)._columnVec;}
        std::vector<double>                     crystalEnergyDepVector(size_t &n){return _vec.at(n)._cryEdepVec;}
        std::vector<double>                     COGrowVector(size_t &n){return _vec.at(n)._COGrowVec;}
        std::vector<double>                     COGcolumnVector(size_t &n){return _vec.at(n)._COGcolumnVec;}
        int                                     vaneId(size_t &n){return _vec.at(n)._vaneId;}
        int                                     COGrow(size_t &n){return _vec.at(n)._cluCogRow;}
        int                                     COGcolumn(size_t &n){return _vec.at(n)._cluCogColumn;}
        int                                     crystalEnergyDepMaxRow(size_t &n){return _vec.at(n)._cryEnergydepMaxRow;}
        int                                     crystalEnergyDepMaxColumn(size_t &n){return _vec.at(n)._cryEnergydepMaxColumn;}
        CLHEP::Hep3Vector                       COG3Vector(size_t &n){return _vec.at(n)._cluCOG;}
        float                                   showerDir(size_t &n){return _vec.at(n)._showerDir;}
        float                                   time(size_t &n){return _vec.at(n)._time;}
        float                                   errShowerDir(size_t &n){return _vec.at(n)._errShowerDir;}
        double                                  clusterEnergy(size_t &n){return _vec.at(n).clusterEnergy();}

        int                                     clusterSize(unsigned int &n){return _vec.at(n)._cluSize;}
        int                                     COGcrystalSize(unsigned int &n){return _vec.at(n)._COGcrySize;}
        std::vector<double>                     rowVector(unsigned int &n){return _vec.at(n)._rowVec;}
        std::vector<double>                     columnVector(unsigned int &n){return _vec.at(n)._columnVec;}
        std::vector<double>                     crystalEnergyDepVector(unsigned int &n){return _vec.at(n)._cryEdepVec;}
        std::vector<double>                     COGrowVector(unsigned int &n){return _vec.at(n)._COGrowVec;}
        std::vector<double>                     COGcolumnVector(unsigned int &n){return _vec.at(n)._COGcolumnVec;}
        int                                     vaneId(unsigned int &n){return _vec.at(n)._vaneId;}
        int                                     COGrow(unsigned int &n){return _vec.at(n)._cluCogRow;}
        int                                     COGcolumn(unsigned int &n){return _vec.at(n)._cluCogColumn;}
        int                                     crystalEnergyDepMaxRow(unsigned int &n){return _vec.at(n)._cryEnergydepMaxRow;}
        int                                     crystalEnergyDepMaxColumn(unsigned int &n){return _vec.at(n)._cryEnergydepMaxColumn;}
        CLHEP::Hep3Vector                       COG3Vector(unsigned int &n){return _vec.at(n)._cluCOG;}
        float                                   showerDir(unsigned int &n){return _vec.at(n)._showerDir;}
        float                                   time(unsigned int &n){return _vec.at(n)._time;}
        float                                   errShowerDir(unsigned int &n){return _vec.at(n)._errShowerDir;}
        double                                  clusterEnergy(unsigned int &n){return _vec.at(n).clusterEnergy();}


        ClusterMap &operator[](unsigned int n){
                return (this->_vec.at(n));
        }
        ClusterMap &at(unsigned int &n){
                return (this->_vec.at(n));
        }

        ClusterMapVector & operator=(const ClusterMapVector other){
                for(size_t i =0; i<other._vec.size(); ++i){
                        _vec.push_back( other._vec.at(i) );
                }
                return *this;
        }

        void print(std::ostream& os) const{
                for(size_t i=0; i<_vec.size(); ++i){
                        os << "ClusterMapVector["<<i<<"] = ";
                        _vec[i].print(os);
                }
        }
        ClusterMapVector(){};

};

//parameters used to build the cluster
class CaloClusterParameters{
public:
        CaloClusterParameters(){}

        double                                _deltaTime;//[ns] time window requested to the crystals of each cluster
        int                             _nCryPerCluster;// minimum number of crystals for defining a cluster
        double                                _EnoiseCut;//[MeV] equal to 3 sigma noise
        double                              _EclusterCut;
        int                                  _nRow;
        int                                _nColum ;

        ~CaloClusterParameters(){}

};

//the following procedure fill the CaloCluster object with the COGVector, and its error vector
// with the respective values obtained using an energy weighted algorithm
void cog(CaloCluster &cluster);

void cog_depth(CaloCluster &cluster, double depth, ClusterMap &clusterMap);

//on the following we implement an algorithm for the cog which uses the logarithm of the energy as weight (w_{i}). This is not correct, as references show, we need to calculate from simulation an offset to add at each w_{i}
//void LOGcog(CaloCluster &cluster);
CLHEP::Hep3Vector LOGcog(CaloCluster &cluster, double w, double depth);

void LOGcogMap(CaloCluster &cluster, double w, double depth, ClusterMap &clusterMap );

 

}

//#endif /* CaloClusterUtilities_hh */
