///////////////////////////////////////////////////////////////////////////////
// TZClusterFilter
// A. M. Ricci and H. Kitagawa
//
// The module works in cartesian and cylindrical coordinates (R, Phi, Z).
///////////////////////////////////////////////////////////////////////////////

// art
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// Offline
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

// ROOT
#include "TH1F.h"

namespace mu2e {

    class TZClusterFilter: public art::EDFilter {

        //-----------------------------------------------------------------------------
        // functions and data structures
        //-----------------------------------------------------------------------------

        public:

        struct Config {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;

            fhicl::Atom<int>           diagLevel          {Name("diagLevel"          ), Comment("turn print on or off"                                 ) };
            fhicl::Atom<int>           runDisplay         {Name("runDisplay"         ), Comment("turn histograms on or off"                            ) };
            fhicl::Atom<art::InputTag> chCollLabel        {Name("chCollLabel"        ), Comment("combo hit collection label"                           ) };
            fhicl::Atom<art::InputTag> tcCollLabel        {Name("tcCollLabel"        ), Comment("time cluster collection Label"                        ) };
            fhicl::Atom<size_t>        minSHsInCluster    {Name("minSHsInCluster"    ), Comment("minimum number of straw hits in the time cluster"     ) };
            fhicl::Atom<size_t>        minCHsInCluster    {Name("minCHsInCluster"    ), Comment("minimum number of combo hits in the time cluster"     ) };
            fhicl::Atom<size_t>        nPairs             {Name("nPairs"             ), Comment("number of pairs of neighboring stations"              ) };
            fhicl::Atom<size_t>        minCHsInCubes      {Name("minCHsInCubes"      ), Comment("minimum number of combo hits in the cubes dR-dPhi-dZ" ) };
            fhicl::Atom<float>         maxDeltaRInStn     {Name("maxDeltaRInStn"     ), Comment("maximum deltaR (mm) in the station"                   ) };
            fhicl::Atom<float>         maxDeltaRBtwStns   {Name("maxDeltaRBtwStns"   ), Comment("maximum deltaR (mm) between stations"                 ) };
            fhicl::Atom<float>         maxDeltaPhiInStn   {Name("maxDeltaPhiInStn"   ), Comment("maximum deltaPhi in the station"                      ) };
            fhicl::Atom<float>         maxDeltaPhiBtwStns {Name("maxDeltaPhiBtwStns" ), Comment("maximum deltaPhi between stations"                    ) };
            fhicl::Atom<bool>          thirdStation       {Name("thirdStation"       ), Comment("check the stations near the pairs"                    ) };
        };

        struct Hist {
            TH1F* nComboHitsPTC;          // Number of combo hits per time cluster
            TH1F* nTClustersPE[2];        // Number of time clusters per event before and after the selection
            TH1F* nTClusters;             // Number of time clusters which pass the selection
        };

        struct comboHit {
            bool center;
            int chCollIdx;
            int nStrawHits;
            int station;
            int plane;
            int face;
            int panel;
            XYZVectorF pos;
        };

        explicit TZClusterFilter(const art::EDFilter::Table<Config>& config);
        virtual ~TZClusterFilter();

        virtual void beginJob();
        virtual bool filter(art::Event& event);

        private:

        //-----------------------------------------------------------------------------
        // helper functions
        //-----------------------------------------------------------------------------
        bool findData(art::Event& event);
        bool fillClusterHits(size_t tc);
        bool checkPopulation(size_t tc);
        bool checkPattern(size_t tc);
        bool checkStation(size_t stn);
        bool checkPair(size_t stn1, size_t stn2);
        bool checkThirdStation(size_t stn1, size_t stn2);
        void clusterInfo(size_t tc);

        //-----------------------------------------------------------------------------
        // functions to manage histograms
        //-----------------------------------------------------------------------------
        void bookHistograms(art::ServiceHandle<art::TFileService>& Tfs);
        void fillComboHitsHistograms(size_t tc);
        void fillHistograms();

        //-----------------------------------------------------------------------------
        // attributes for diagnostics
        //-----------------------------------------------------------------------------
        int                 _diagLevel;
        int                 _runDisplay;
        Hist                _hist;
        art::Event*         _event;
        std::vector<size_t> _nTCs;

        //-----------------------------------------------------------------------------
        // attributes for event object labels
        //-----------------------------------------------------------------------------
        art::InputTag _chLabel;
        art::InputTag _tcLabel;

        //-----------------------------------------------------------------------------
        // attributes for collections
        //-----------------------------------------------------------------------------
        const ComboHitCollection*    _chColl;
        const TimeClusterCollection* _tcColl;
        TimeClusterCollection*       _tcColl2;

        //-----------------------------------------------------------------------------
        // attributes for cluster selection parameters
        //-----------------------------------------------------------------------------
        size_t _minSHsInCluster;
        size_t _minCHsInCluster;
        size_t _nPairs;
        size_t _minCHs;
        float  _deltaR1;
        float  _deltaR2;
        float  _deltaPhi_1;
        float  _deltaPhi_2;
        bool   _thirdStation;

        //-----------------------------------------------------------------------------
        // attributes used for cluster selection
        //-----------------------------------------------------------------------------
        std::vector<comboHit> _tcHits;
        std::vector<comboHit> _chStn[StrawId::_nstations];
    };

    //-----------------------------------------------------------------------------
    // module constructor
    //-----------------------------------------------------------------------------
    TZClusterFilter::TZClusterFilter(const art::EDFilter::Table<Config>& config) :
    art::EDFilter(config),
    _diagLevel              (config().diagLevel()                                 ),
    _runDisplay             (config().runDisplay()                                ),
    _chLabel                (config().chCollLabel()                               ),
    _tcLabel                (config().tcCollLabel()                               ),
    _minSHsInCluster        (config().minSHsInCluster()                           ),
    _minCHsInCluster        (config().minCHsInCluster()                           ),
    _nPairs                 (config().nPairs()                                    ),
    _minCHs                 (config().minCHsInCubes()                             ),
    _deltaR1                (config().maxDeltaRInStn()                            ),
    _deltaR2                (config().maxDeltaRBtwStns()                          ),
    _deltaPhi_1             (config().maxDeltaPhiInStn()                          ),
    _deltaPhi_2             (config().maxDeltaPhiBtwStns()                        ),
    _thirdStation           (config().thirdStation()                              )
    {
        // tell art the data products that the module will use and produce
        consumes<ComboHitCollection>(_chLabel);
        consumes<TimeClusterCollection>(_tcLabel);
        produces<TimeClusterCollection>();
    }

    //-----------------------------------------------------------------------------
    // destructor
    //-----------------------------------------------------------------------------
    TZClusterFilter::~TZClusterFilter() {}

    //-----------------------------------------------------------------------------
    // beginJob
    //-----------------------------------------------------------------------------
    void TZClusterFilter::beginJob() {
        if(_runDisplay > 0) {
            art::ServiceHandle<art::TFileService> tfs;
            bookHistograms(tfs);
        }
    }

    //-----------------------------------------------------------------------------
    // book histograms
    //-----------------------------------------------------------------------------
    void TZClusterFilter::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

        TH1::AddDirectory(0);

        _hist.nComboHitsPTC = Tfs->make<TH1F>("nComboHitsPTC0", "Number of combo hits per time cluster", 150, 0.0, 150.0 );
        _hist.nTClustersPE[0] = Tfs->make<TH1F>("nTClustersPE0", "Number of time clusters per event before the selection", 31, 0.0, 31.0 );
        _hist.nTClustersPE[1] = Tfs->make<TH1F>("nTClustersPE1", "Number of time clusters per event after the selection", 31, 0.0, 31.0 );
        _hist.nTClusters = Tfs->make<TH1F>("nTClusters", "Number of time clusters which pass the selection", 3, 0.0, 3.0 );
    }

    //-----------------------------------------------------------------------------
    // fill histograms containing information on combo hits
    //-----------------------------------------------------------------------------
    void TZClusterFilter::fillComboHitsHistograms(size_t tc) {

        // fill the number of combo hits per time cluster
        _hist.nComboHitsPTC->Fill(_tcHits.size());
    }

    //-----------------------------------------------------------------------------
    // fill histograms containing information on time clusters
    //-----------------------------------------------------------------------------
    void TZClusterFilter::fillHistograms() {

        // fill the number of time clusters before the selection
        _hist.nTClustersPE[0]->Fill(_tcColl->size());

        // fill the number of time clusters after the selection
        _hist.nTClustersPE[1]->Fill(_tcColl2->size());

        // fill the number of time clusters which pass the selection
        for(size_t i=0; i<_nTCs.size(); i++) _hist.nTClusters->Fill(_nTCs[i]);
    }

    //-----------------------------------------------------------------------------
    // find input
    //-----------------------------------------------------------------------------
    bool TZClusterFilter::findData(art::Event& event) {

        // get Combo Hit Collection
        auto chcolH = event.getValidHandle<ComboHitCollection>(_chLabel);
        if(chcolH.product() != 0) {
            _chColl = chcolH.product();
        }
        else {
            _chColl = 0;
            std::cout << ">>> ERROR in TZClusterFilter::findData: ComboHitCollection not found." << std::endl;
        }

        // get Time Cluster Collection
        auto tccolH = event.getValidHandle<TimeClusterCollection>(_tcLabel);
        if(tccolH.product() != 0) {
            _tcColl = tccolH.product();
        }
        else {
            _tcColl = 0;
            std::cout << ">>> ERROR in TZClusterFilter::findData: TimeClusterCollection not found." << std::endl;
        }

        return (_chColl != 0 && _tcColl != 0);
    }

    //-----------------------------------------------------------------------------
    // fill Combo Hits of the Time Clusters
    //-----------------------------------------------------------------------------
    bool TZClusterFilter::fillClusterHits(size_t tc) {

        _tcHits.clear();

        comboHit cHit;

        // loop on the combo hits of the time cluster
        // _strawHitIdxs is the index of the combo hits in the time cluster
        for(size_t i=0; i<_tcColl->at(tc)._strawHitIdxs.size(); i++) {
            int hitIndex = _tcColl->at(tc)._strawHitIdxs[i];
            cHit.center = false;
            cHit.chCollIdx = hitIndex;
            cHit.nStrawHits = _chColl->at(hitIndex).nStrawHits();
            cHit.station = _chColl->at(hitIndex).strawId().station();
            cHit.plane = _chColl->at(hitIndex).strawId().plane();
            cHit.face = _chColl->at(hitIndex).strawId().face();
            cHit.panel = _chColl->at(hitIndex).strawId().panel();
            cHit.pos = _chColl->at(hitIndex).pos();

            _tcHits.push_back(cHit);
        }

        // sort the vector in ascending order of the z-coordinate or print an error message
        if(_tcHits.empty()) std::cout << ">>> ERROR in TZClusterFilter::fillTClusterHits: the vector _tcHits is empty." << std::endl;
        else std::sort(_tcHits.begin(), _tcHits.end(), [](const comboHit& a, const comboHit& b) { return a.pos.z() < b.pos.z(); } );

        return (_tcHits.empty());
    }

    //-----------------------------------------------------------------------------
    // event entry point
    //-----------------------------------------------------------------------------
    bool TZClusterFilter::filter(art::Event& event) {

        _nTCs.clear();

        _event = &event;
        bool data = findData(event);

        if(!data) return false;

        std::unique_ptr<TimeClusterCollection> tcColl(new TimeClusterCollection);
        _tcColl2 = tcColl.get();

        // loop on Time Clusters
        bool cluster = false;
        for(size_t i=0; i<_tcColl->size(); i++) {
            bool tcHits = fillClusterHits(i);
            if(tcHits) {
                _nTCs.push_back(0);
                continue;
            }
            if(_diagLevel > 0) clusterInfo(i);
            if(_runDisplay > 0) fillComboHitsHistograms(i);
            bool popCheck = checkPopulation(i);
            if(!popCheck) {
                _nTCs.push_back(0);
                continue;
            }
            bool pattern = checkPattern(i);
            if(!pattern) {
                _nTCs.push_back(0);
                continue;
            }
            _tcColl2->push_back(_tcColl->at(i));
            _nTCs.push_back(1);
            cluster = true;
        }

        if(_runDisplay > 0) fillHistograms();

        //-----------------------------------------------------------------------------
        // put time cluster collection into the event record
        //-----------------------------------------------------------------------------
        event.put(std::move(tcColl));

        return cluster;
    }

    //-----------------------------------------------------------------------------
    // check the population of the time cluster
    //-----------------------------------------------------------------------------
    bool TZClusterFilter::checkPopulation(size_t tc) {

        size_t TotalStrawHits = 0;
        for(size_t i=0; i<_tcHits.size(); i++) {
            TotalStrawHits += _tcHits[i].nStrawHits;
        }

        // the time cluster is rejected if it has too few straw hits
        if(TotalStrawHits < _minSHsInCluster) {
            if(_diagLevel > 0) {
                std::cout << ">>> WARNING in TZClusterFilter::checkPopulation:" << std::endl;
                std::cout << ">>> The Time Cluster " << tc << " contains less than " << _minSHsInCluster <<
                                                                " Straw Hits and is rejected." << std::endl;
            }
            return false;
        }

        // the time cluster is rejected if it has too few combo hits
        if(_tcHits.size() < _minCHsInCluster) {
            if(_diagLevel > 0) {
                std::cout << ">>> WARNING in TZClusterFilter::checkPopulation:" << std::endl;
                std::cout << ">>> The Time Cluster " << tc << " contains less than " << _minCHsInCluster <<
                                                                " Combo Hits and is rejected." << std::endl;
            }
            return false;
        }

        return true;
    }

    //-----------------------------------------------------------------------------
    // recognize the pattern in the time cluster
    //-----------------------------------------------------------------------------
    bool TZClusterFilter::checkPattern(size_t tc) {

        for(size_t i=0; i<StrawId::_nstations; i++) _chStn[i].clear();

        // divide the combo hits per station
        for(size_t i=0; i<_tcHits.size(); i++) {
            for(size_t j=0; j<StrawId::_nstations; j++) {
                if(_tcHits[i].station == int(j)) _chStn[j].push_back(_tcHits[i]);
            }
        }

        // loop on the stations
        std::vector<size_t> pairs;
        for(size_t i=0; i<17; i++) {
            // take the pairs of neighboring stations with at least _minCHs combo hits per station
            if(_chStn[i].size() < _minCHs || _chStn[i+1].size() < _minCHs) continue;

            if(_diagLevel > 1) {
                std::cout << ">>> INFORMATION in TZClusterFilter::checkPattern:" << std::endl;
                std::cout << ">>> Candidate pair: Time Cluster " << tc << ", Stations " << i << "/" << i+1 <<
                        " contain " <<_chStn[i].size() << "/" << _chStn[i+1].size() << " Combo Hits." << std::endl;
            }

            bool PairCheck[4] = {false};

            // check whether each station of the pair has at least one
            // rectangle of side _deltaPhi_1 that contains at least 2 CHs
            PairCheck[0] = checkStation(i);
            if(!PairCheck[0]) continue;
            PairCheck[1] = checkStation(i+1);
            if(!PairCheck[1]) continue;

            // check whether the pair of neighboring stations has at least
            // two rectangles within _deltaPhi_2 (one rectangle per station)
            PairCheck[2] = checkPair(i, i+1);
            if(!PairCheck[2]) continue;

            pairs.push_back(i);
        }

        // loop on the pairs of neighboring stations: check whether near a pair
        // there is a third station with at least a combo hit within _deltaPhi_2
        bool TSCheck = false;
        if(_thirdStation && !pairs.empty()) {
            for(size_t i=0; i<pairs.size(); i++) {
                TSCheck = checkThirdStation(pairs[i], pairs[i]+1);
                if(TSCheck) break;
            }
        }

        if(_diagLevel > 0) {
            std::cout << ">>> INFORMATION in TZClusterFilter::checkPattern:" << std::endl;
            std::cout << ">>> Found " << pairs.size() << " pairs ( ";
            for(size_t i=0; i<pairs.size(); i++) std::cout << pairs[i] << "-" << pairs[i]+1 << " ";
            std::cout << ") in Time Cluster " << tc << "." << std::endl;
        }

        if(pairs.size() < _nPairs) return false;
        if(_thirdStation && !TSCheck) return false;

        return true;
    }

    //-----------------------------------------------------------------------------
    // check the combo hits in the station
    //-----------------------------------------------------------------------------
    bool TZClusterFilter::checkStation(size_t stn) {

        bool rectangle = false;
        size_t nCHs = 0;

        // combinatorial computation: for each combo hit on the station count
        // the number of combo hits inside the cube of sides: _deltaR1 and _deltaPhi_1
        for(size_t i=0; i<_chStn[stn].size(); i++) {
            for(size_t j=0; j<_chStn[stn].size(); j++) {
                // same hit
                if(i==j) {
                    nCHs++;
                    continue;
                }
                double a[2] = {_chStn[stn].at(i).pos.x(), _chStn[stn].at(i).pos.y()};
                double b[2] = {_chStn[stn].at(j).pos.x(), _chStn[stn].at(j).pos.y()};
                double Ra = std::sqrt(a[0]*a[0]+a[1]*a[1]);
                double Rb = std::sqrt(b[0]*b[0]+b[1]*b[1]);
                double deltaR = std::abs(Rb - Ra);
                if(deltaR > _deltaR1) continue;
                // scalar product
                double c = (a[0]*b[0]+a[1]*b[1])/(Ra*Rb);
                double deltaPhi = acos(c); // [rad]
                if(deltaPhi <= _deltaPhi_1) nCHs++;
            }
            if(_diagLevel > 1) {
                std::cout << ">>> INFORMATION in TZClusterFilter::checkStation:" << std::endl;
                std::cout << ">>> Station " << stn << ", Cubes: CentralCH(chCollIdx)/nCHs = " <<
                                                i << "(" << _chStn[stn].at(i).chCollIdx << ")/" << nCHs << std::endl;
            }
            if(nCHs >= _minCHs) {
                _chStn[stn].at(i).center = true;
                rectangle = true;
            }
            nCHs = 0;
        }

        return rectangle;
    }

    //-----------------------------------------------------------------------------
    // check the pair of neighboring stations
    //-----------------------------------------------------------------------------
    bool TZClusterFilter::checkPair(size_t stn1, size_t stn2) {

        // combinatorial computation: check whether the pair of neighboring stations
        // has at least two cubes within _deltaR2 and _deltaPhi_2 (one cube per station)
        // note: the centers of the cubes are combo hits.
        for(size_t i=0; i<_chStn[stn1].size(); i++) {
            if(_chStn[stn1].at(i).center == false) continue;
            for(size_t j=0; j<_chStn[stn2].size(); j++) {
                if(_chStn[stn2].at(j).center == false) continue;
                double a[2] = {_chStn[stn1].at(i).pos.x(), _chStn[stn1].at(i).pos.y()};
                double b[2] = {_chStn[stn2].at(j).pos.x(), _chStn[stn2].at(j).pos.y()};
                double Ra = std::sqrt(a[0]*a[0]+a[1]*a[1]);
                double Rb = std::sqrt(b[0]*b[0]+b[1]*b[1]);
                double deltaR = std::abs(Rb - Ra);
                if(deltaR > _deltaR2) continue;
                // scalar product
                double c = (a[0]*b[0]+a[1]*b[1])/(Ra*Rb);
                double deltaPhi = acos(c); // [rad]
                if(deltaPhi <= _deltaPhi_2) {
                    if(_diagLevel > 1) {
                        std::cout << ">>> INFORMATION in TZClusterFilter::checkPair:" << std::endl;
                        std::cout << ">>> Found pair: Stations " << stn1 << "/" << stn2 <<
                                        ", Cubes: CentralCH(chCollIdx)/CentralCH(chCollIdx) = " <<
                                                    i << "(" << _chStn[stn1].at(i).chCollIdx << ")/" <<
                                                    j << "(" << _chStn[stn2].at(j).chCollIdx << ")" << std::endl;
                    }
                    return true;
                }
            }
        }

        return false;
    }

    //-----------------------------------------------------------------------------
    // check the third station near the pair of neighboring stations
    //-----------------------------------------------------------------------------
    bool TZClusterFilter::checkThirdStation(size_t stn1, size_t stn2) {

        // combinatorial computation: check whether near the pair of neighboring stations
        // there is a third station with at least a combo hit within _deltaR2 and _deltaPhi_2
        // note: the centers of the cubes are combo hits.
        if(stn1 > 0) {
            for(size_t i=0; i<_chStn[stn1].size(); i++) {
                if(_chStn[stn1].at(i).center == false) continue;
                for(size_t j=0; j<_chStn[stn1-1].size(); j++) {
                    double a[2] = {_chStn[stn1].at(i).pos.x(), _chStn[stn1].at(i).pos.y()};
                    double b[2] = {_chStn[stn1-1].at(j).pos.x(), _chStn[stn1-1].at(j).pos.y()};
                    double Ra = std::sqrt(a[0]*a[0]+a[1]*a[1]);
                    double Rb = std::sqrt(b[0]*b[0]+b[1]*b[1]);
                    double deltaR = std::abs(Rb - Ra);
                    if(deltaR > _deltaR2) continue;
                    // scalar product
                    double c = (a[0]*b[0]+a[1]*b[1])/(Ra*Rb);
                    double deltaPhi = acos(c); // [rad]
                    if(deltaPhi <= _deltaPhi_2) {
                        if(_diagLevel > 0) {
                            std::cout << ">>> INFORMATION in TZClusterFilter::checkThirdStation:" << std::endl;
                            std::cout << ">>> Found third station before the pair " << stn1 << "/" << stn2 <<
                                ", ComboHit/chCollIdx: " << j << "/" << _chStn[stn1-1].at(j).chCollIdx << std::endl;
                        }
                        return true;
                    }
                }
            }
        }

        if(stn2 < 17) {
            for(size_t i=0; i<_chStn[stn2].size(); i++) {
                if(_chStn[stn2].at(i).center == false) continue;
                for(size_t j=0; j<_chStn[stn2+1].size(); j++) {
                    double a[2] = {_chStn[stn2].at(i).pos.x(), _chStn[stn2].at(i).pos.y()};
                    double b[2] = {_chStn[stn2+1].at(j).pos.x(), _chStn[stn2+1].at(j).pos.y()};
                    double Ra = std::sqrt(a[0]*a[0]+a[1]*a[1]);
                    double Rb = std::sqrt(b[0]*b[0]+b[1]*b[1]);
                    double deltaR = std::abs(Rb - Ra);
                    if(deltaR > _deltaR2) continue;
                    // scalar product
                    double c = (a[0]*b[0]+a[1]*b[1])/(Ra*Rb);
                    double deltaPhi = acos(c); // [rad]
                    if(deltaPhi <= _deltaPhi_2) {
                        if(_diagLevel > 0) {
                            std::cout << ">>> INFORMATION in TZClusterFilter::checkThirdStation:" << std::endl;
                            std::cout << ">>> Found third station after the pair " << stn1 << "/" << stn2 <<
                                ", ComboHit/chCollIdx: " << j << "/" << _chStn[stn2+1].at(j).chCollIdx << std::endl;
                        }
                        return true;
                    }
                }
            }
        }

        return false;
    }

    //-----------------------------------------------------------------------------
    // print information on the time cluster
    //-----------------------------------------------------------------------------
    void TZClusterFilter::clusterInfo(size_t tc) {

        std::cout << ">>> INFORMATION in TZClusterFilter::clusterInfo: " << std::endl;
        std::cout << "===========================================" << std::endl;
        std::cout << " Time Cluster " << tc << std::endl;
        std::cout << " Total Combo Hits: " << _tcColl->at(tc)._strawHitIdxs.size() << std::endl;
        std::cout << " Combo Hits selected: " << _tcHits.size() << std::endl;
        std::cout << "===========================================" << std::endl;

        if(_diagLevel < 2) return;

        for(size_t i=0; i<_tcHits.size(); i++) {
            std::cout << "ComboHit/chCollIdx: " << i << "/" << _tcHits[i].chCollIdx << std::endl;
            std::cout << "nStrawHits/station/plane/face/panel: " <<
                _tcHits[i].nStrawHits << "/" << _tcHits[i].station << "/" << _tcHits[i].plane <<
                                        "/" << _tcHits[i].face << "/" << _tcHits[i].panel << std::endl;
            std::cout << "x/y/z/phi: " << _tcHits[i].pos.x() << "/" << _tcHits[i].pos.y() <<
                            "/" << _tcHits[i].pos.z() << "/" << _tcHits[i].pos.phi() << std::endl;
        }
    }

}

using mu2e::TZClusterFilter;
DEFINE_ART_MODULE(TZClusterFilter)
