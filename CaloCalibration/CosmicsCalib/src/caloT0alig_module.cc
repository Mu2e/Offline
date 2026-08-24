// =============================================================================
// Filter module for offline calorimeter time calibration
// 18-Mar-2021 - S.Giovannella
//
// Starting code for retrieving calo information: CaloExample module
//
// This module performs a single step of an iterative procedure running on
// calorimeter readout channels:
// 1. selection of straight cosmic ray events with mip-like energy deposit
// 2. fit to the 2d linear trajectory in the XY plane with least-square method
// 3. common energy-weighted T0 subtracted to all readout channel
// 4. speed of light imposed between timing of different readouts
// 5. residuals evaluated for each readout channel, to be used as input for
//    the next iteration
// Iterations are controlled by an external script
// Good events can be filtered to speed up the procedure
// Last iteration creates an output file to be uploaded on Condition DB
// =============================================================================

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"

#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"

#include <fstream>
#include <iostream>
#include <string>

#include "TDirectory.h"
#include "TH1F.h"
#include <TF1.h>
#include <TGraphErrors.h>

namespace mu2e {

class caloT0alig : public art::EDFilter {

public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<art::InputTag> caloHitCollection{Name("caloHitCollection"),
                                                 Comment("Calo Hit collection name")};
    fhicl::Atom<art::InputTag> caloClusterCollection{Name("caloClusterCollection"),
                                                     Comment("Calo cluster collection name")};
    fhicl::Atom<int> cluHits{Name("cluHits"), Comment("Minimum number of hits in a cluster"), 3};
    fhicl::Atom<double> cryEmin{Name("cryEmin"), Comment("Minimum Energy for crystal hits"), 15.};
    fhicl::Atom<int> ncryCut{Name("ncryCut"), Comment("Minimum number of mip-like crystals"), 7};
    fhicl::Atom<std::string> iteration{Name("iteration"), Comment("Iteration"), "first"};
    fhicl::Atom<std::string> fileT0{Name("fileT0"), Comment("T0 input file")};
    fhicl::Atom<std::string> fileTcor{Name("fileTcor"), Comment("T0 corrections input file")};
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("Diag Level"), 0};
  };

  explicit caloT0alig(const art::EDFilter::Table<Config>& config);
  virtual ~caloT0alig() {}

  virtual void beginJob();
  virtual void endJob();
  virtual bool filter(art::Event& e);

private:
  art::InputTag _caloHitTag;
  art::InputTag _caloClusterTag;
  int _cluHits;
  float _cryEmin;
  int _ncryCut;
  std::string _iteration;
  std::ifstream _fileT0;
  std::string _fileT0Name;
  std::string _fileTcorName;
  int _diagLevel;
  int _nProcessed;
  int _nFiltered;

  static constexpr double cvel = 299.792458;

  static constexpr int nDisks = CaloConst::_nDisk;              // Number of calorimeter disks
  static constexpr int nCrystals = CaloConst::_nCrystalPerDisk; // Number of crystals for each disk
  static constexpr int nSiPMs = CaloConst::_nSiPMPerCrystal;    // Number of SIPMs for each crystal
  static constexpr int nROchan = nDisks * nCrystals * nSiPMs;

  float Tcor[nROchan] = {0};
  float Toff[nROchan] = {0};

  TH1F* hTres[nROchan];
  TH1F *hcosTh[nDisks], *hcosThCut[nDisks];
  TH1F *hcosFit[nDisks], *hChi2[nDisks], *hDcosth[nDisks], *hYdif[nDisks];
  TH1F *hcvel[nDisks], *hEcell[nDisks], *hEmean[nDisks];
  TH1F *hNphi[nDisks], *hNphiEcut[nDisks];
  TH1F *hNhit, *hNhitEcut;
};

caloT0alig::caloT0alig(const art::EDFilter::Table<Config>& config) :
    EDFilter{config}, _caloHitTag(config().caloHitCollection()),
    _caloClusterTag(config().caloClusterCollection()), _cluHits(config().cluHits()),
    _cryEmin(config().cryEmin()), _ncryCut(config().ncryCut()), _iteration(config().iteration()),
    _fileT0Name(config().fileT0()), _fileTcorName(config().fileTcor()), _diagLevel(config().diagLevel()),
    _nProcessed(0), _nFiltered(0) {

  _fileT0(_fileT0Name());
  if (!_fileT0.is_open()) {
    throw cet::exception("caloT0alig")
        << "ERROR! Cannot open input file " << config().fileT0() << std::endl;
  }
}

// ===========================================================================
// Begin job:
// - Get first step of T0 corrections, i.e. mean values of time distribution
//   for cosmic ray events
// - Get T0 corrections from previous iterations, if applicable
// - Histogram booking
// ===========================================================================
void caloT0alig::beginJob() {

  if (_diagLevel > 0)
    std::cout << "caloT0alig: Entering beginJob" << std::endl;

  // Read T0-step0 corrections
  // *** TO BE IMPLEMENTED ... waiting for real data ***

  ///////////////////////////////////////////////////////////////////////
  // TEMPORARY for running on MC data:
  // T0 smearing (+-0.5 ns) to simulate real data and test the procedure
  ///////////////////////////////////////////////////////////////////////

  int iChanT0, nValT0 = 0;
  float TvalT0;

  if (_fileT0.is_open()) {
    while (_fileT0 >> iChanT0 >> TvalT0) {
      if (iChanT0 < 0 || iChanT0 >= nROchan) {
        throw cet::exception("caloT0alig") << "ERROR! Read invalid channel " << iChanT0
                                           << "from file " << _fileT0Name << std::endl;
      }
      Toff[iChanT0] = TvalT0;
      if (_diagLevel > 1)
        std::cout << "IdxT0 " << iChanT0 << " Toff " << Toff[iChanT0] << " " << std::endl;
      nValT0++;
    }
    _fileT0.close();
  } else {
    mf::LogError("INPUT-NOT-FOUND")
        << "T0 file from previous iteration not found: " << _fileT0Name << std::endl;
  }
  ///////////////////////////////////////////////////////////////

  // Read T0 corrections from previous iteration, if applicable

  if (_iteration == "middle" || _iteration == "last") {
    if (_diagLevel > 0)
      std::cout << "Iteration " << _iteration << ": Reading T0 corrections from previous step "
                << std::endl;

    int iChan, nVal = 0;
    int Nev;
    float Tval, Tmea, Tres, Chi2;

    std::ifstream _fileTcor(_fileTcorName);
    if (_fileTcor.is_open()) {
      while (_fileTcor >> iChan >> Tval >> Tmea >> Tres >> Chi2 >> Nev) {
        if (iChan < 0 || iChan >= nROchan) {
          throw cet::exception("caloT0alig") << "ERROR! Read invalid channel " << iChan
                                             << "from file " << _fileTcorName << std::endl;
        }

        Tcor[iChan] = Tval;
        if (_diagLevel > 1)
          std::cout << "Idx " << iChan << " Tcor " << Tcor[iChan] << " " << Tmea << std::endl;
        nVal++;
      }
      _fileTcor.close();
      if (nVal != nROchan) {
        mf::LogError("WRONG-NCHAN")
            << "Wrong number of readout channels: " << nVal << " " << nROchan << std::endl;
      }
    } else {
      if (!_fileTcor.is_open()) {
        throw cet::exception("caloT0alig")
            << "ERROR! Cannot open output file " << _fileTcorName << std::endl;
      }
    }
  }

  // Book histograms

  art::ServiceHandle<art::TFileService> tfs;

  for (int iChan = 0; iChan < nROchan; iChan++) {
    hTres[iChan] =
        tfs->make<TH1F>(Form("Tres_%04i", iChan), Form("Tres %04i", iChan), 200, -2., 2.);
  }

  for (int iDisk = 0; iDisk < nDisks; iDisk++) {
    hcosTh[iDisk] = tfs->make<TH1F>(Form("cosTh_%i", iDisk), Form("cosTh %i", iDisk), 200, -1., 1.);
    hcosThCut[iDisk] =
        tfs->make<TH1F>(Form("cosThCut_%i", iDisk), Form("cosThCut %i", iDisk), 200, -1., 1.);
    hcosFit[iDisk] =
        tfs->make<TH1F>(Form("cosFit_%i", iDisk), Form("cosFit %i", iDisk), 200, -1., 1.);
    hDcosth[iDisk] =
        tfs->make<TH1F>(Form("Dcosth_%i", iDisk), Form("Dcosth %i", iDisk), 2000, -1., 1.);
    hYdif[iDisk] =
        tfs->make<TH1F>(Form("Ydif_%i", iDisk), Form("Ydif %i", iDisk), 200, -1000., 1000.);
    hChi2[iDisk] = tfs->make<TH1F>(Form("Chi2_%i", iDisk), Form("Chi2 %i", iDisk), 200, 0., 100.);
    hcvel[iDisk] = tfs->make<TH1F>(Form("cvel_%i", iDisk), Form("cvel %i", iDisk), 200, 100., 500.);
    hEcell[iDisk] =
        tfs->make<TH1F>(Form("Ecell_%01i", iDisk), Form("Ecell %i", iDisk), 200, 0., 100.);
    hEmean[iDisk] =
        tfs->make<TH1F>(Form("Emean_%01i", iDisk), Form("Emean %i", iDisk), 200, 0., 100.);

    hNphi[iDisk] =
        tfs->make<TH1F>(Form("Nphi_%01i", iDisk), Form("Nphi %i", iDisk), 180, -180., 180.);
    hNphiEcut[iDisk] =
        tfs->make<TH1F>(Form("NphiEcut_%01i", iDisk), Form("NphiEcut %i", iDisk), 180, -180., 180.);
  }
  hNhit = tfs->make<TH1F>("Nhit", "Nhit", 2 * nCrystals, 0., 2 * nCrystals);
  hNhitEcut = tfs->make<TH1F>("NhitEcut", "NhitEcut", 2 * nCrystals, 0., 2 * nCrystals);

} // End of beginJob

// ===========================================================================
// Event:
// - Select straight cosmic ray tracks crossing calo disks and fill time
//   histogram for each readout channel
// - Filters good events
// ===========================================================================
bool caloT0alig::filter(art::Event& event) {

  bool retval(false); // preset to fail

  ++_nProcessed;

  if (_nProcessed % 10000 == 0 && _diagLevel > 0)
    std::cout << "caloT0alig: Processing run/event " << event.run() << " " << event.id().event()
              << std::endl;

  //
  // *** TO BE IMPLEMENTED ***
  // Check trigger bits to run only on trigger selected cosmics
  //

  // Handle to the calorimeter
  art::ServiceHandle<GeometryService> geom;
  if (!geom->hasElement<Calorimeter>()) {
    mf::LogError("NO-CALO-GEOM") << "Calorimeter geometry not found!" << std::endl;
  }
  const Calorimeter& caloGeom = *(GeomHandle<Calorimeter>());

  // Calorimeter crystal hits (average from readouts)
  art::Handle<CaloHitCollection> CaloHitsHandle;
  event.getByLabel(_caloHitTag, CaloHitsHandle);
  const CaloHitCollection& CaloHits(*CaloHitsHandle);

  // Calorimeter clusters
  art::Handle<CaloClusterCollection> caloClustersHandle;
  event.getByLabel(_caloClusterTag, caloClustersHandle);
  const CaloClusterCollection& caloClusters(*caloClustersHandle);

  float sx[nDisks] = {};
  float sy[nDisks] = {};
  float sxy[nDisks] = {};
  float sy2[nDisks] = {};

  int Ival[nDisks][nROchan];
  float Xval[nDisks][nROchan];
  float Yval[nDisks][nROchan];
  float Eval[nDisks][nROchan];
  float Tval[nDisks][nROchan];
  int nCry[nDisks] = {};  // Good cells on calo disks
  int nChan[nDisks] = {}; // Good readout channels on calo disks

  float discr, bb, costh, Tres;

  for (unsigned int iClu = 0; iClu < caloClusters.size(); ++iClu) { // Loop on clusters
    const CaloCluster& cluster = caloClusters.at(iClu);

    if (cluster.size() >= _cluHits) { // Hits in a cluster
      if (_diagLevel > 2)
        std::cout << "Nproc cluId " << _nProcessed << " " << iClu << std::endl;
      int diskId = cluster.diskID();

      std::vector<int> cryList;
      for (auto cryPtr : cluster.caloHitsPtrVector())
        cryList.push_back(std::distance(&CaloHits.at(0), cryPtr.get()));

      for (int iCry = 0; iCry < cluster.size(); iCry++) { // Crystals connected to cluster

        const CaloHit& hit = CaloHits.at(cryList[iCry]);

        int cryId = hit.crystalID();
        int cryDisk = caloGeom.crystal(cryId).diskID();
        if (cryDisk != diskId)
          mf::LogError("WRONG-DISK")
              << "Different clu/cry disk: " << diskId << " " << cryDisk << std::endl;
        float cryE = hit.energyDep();

        // Hit distribution vs phi and cristal index before t0alig selection cuts

        CLHEP::Hep3Vector crystalPos = caloGeom.mu2eToDiskFF(
            diskId, caloGeom.crystal(cryId).position()); // Crystal position in disk FF frame
        float xval = crystalPos.x();
        float yval = crystalPos.y();

        float phiDeg = atan2(yval, xval) * 180. / 3.1416;
        hNphi[cryDisk]->Fill(phiDeg);
        hNhit->Fill(cryId);
        if (cryE > _cryEmin) {
          hNphiEcut[cryDisk]->Fill(phiDeg);
          hNhitEcut->Fill(cryId);
        }

        if (cryDisk == diskId && cryE > _cryEmin && cryE < 30.) {

          sx[diskId] = sx[diskId] + xval;
          sy[diskId] = sy[diskId] + yval;
          sxy[diskId] = sxy[diskId] + xval * yval;
          sy2[diskId] = sy2[diskId] + pow(yval, 2);

          // Get values for the two SiPMs from digis
          for (unsigned int iCha = 0; iCha < hit.recoCaloDigis().size(); iCha++) {

            int idx = CaloSiPMId(hit.recoCaloDigis().at(iCha)->SiPMID()).SiPMLocalId();
            Ival[diskId][nChan[diskId]] = hit.recoCaloDigis().at(idx)->SiPMID();
            Xval[diskId][nChan[diskId]] = xval;
            Yval[diskId][nChan[diskId]] = yval;
            Eval[diskId][nChan[diskId]] = cryE;
            Tval[diskId][nChan[diskId]] = hit.recoCaloDigis().at(idx)->time();
            Tval[diskId][nChan[diskId]] =
                Tval[diskId][nChan[diskId]] + Toff[Ival[diskId][nChan[diskId]]];
            if (_diagLevel > 2)
              std::cout << "Disk Cry x y idx E T : " << cryDisk << " " << cryId << " " << xval
                        << " " << yval << " " << Ival[diskId][nChan[diskId]] << " " << cryE << " "
                        << Tval[diskId][nChan[diskId]] << std::endl;
            nChan[diskId]++;
          }
          nCry[diskId]++;
        } // Disk + energy cuts
      } // Loop on crystals
    } // Minimum number of cells in a cluster
  } // Loop on clusters

  for (int iDisk = 0; iDisk < nDisks; iDisk++) {

    int Nfit = 0;
    float enetot = 0.;
    double t0 = 0.;
    float Xfit[nCrystals] = {}, Yfit[nCrystals] = {}, Efit[nCrystals] = {};

    if (nCry[iDisk] > _ncryCut) { // Enough good cells

      discr = float(nCry[iDisk]) * sy2[iDisk] - pow(sy[iDisk], 2);

      if (discr != 0.) {

        bb = (float(nCry[iDisk]) * sxy[iDisk] - sx[iDisk] * sy[iDisk]) / discr;
        costh = 1. / sqrt(1. + bb * bb);
        hcosTh[iDisk]->Fill(costh);

        if (costh > 0.2) {
          hcosThCut[iDisk]->Fill(costh);

          for (int iCha = 0; iCha < nChan[iDisk]; iCha++) {
            t0 = t0 + Eval[iDisk][iCha] * (Tval[iDisk][iCha] - Tcor[Ival[iDisk][iCha]] +
                                           Yval[iDisk][iCha] / (cvel * costh));
            enetot = enetot + Eval[iDisk][iCha];
            hEcell[iDisk]->Fill(Eval[iDisk][iCha]);

            Xfit[Nfit] = Xval[iDisk][iCha];
            Yfit[Nfit] = Yval[iDisk][iCha];
            Efit[Nfit] = 10.;
            if (_diagLevel > 3)
              std::cout << Nfit << " " << Xfit[Nfit] << " " << Yfit[Nfit] << " " << Efit[Nfit]
                        << std::endl;
            Nfit++;
          }

          // Linear fit to x,y cell positions
          TGraphErrors* emcpos = new TGraphErrors(Nfit, Xfit, Yfit, Efit, Efit);
          TF1* linfit = new TF1("linfit", "[0]+[1]*x", -640., 640.);
          emcpos->Fit("linfit", "rboq");
          double Par0 = linfit->GetParameter(0);
          double Par1 = linfit->GetParameter(1);
          double Chi2 = linfit->GetChisquare();
          if (_diagLevel > 3)
            std::cout << Par0 << " " << Par1 << std::endl;
          if (_diagLevel > 3)
            std::cout << "chi2: " << Chi2 << std::endl;
          Chi2 = Chi2 / linfit->GetNDF();
          hChi2[iDisk]->Fill(Chi2);
          float costhFit = cos(atan(Par1));
          hcosFit[iDisk]->Fill(costhFit);

          if (Chi2 < 2. && enetot != 0.) {

            t0 = t0 / enetot;
            hEmean[iDisk]->Fill(enetot / nChan[iDisk]);

            // Good event
            retval = true;
            ++_nFiltered;

            for (int iCha = 0; iCha < nChan[iDisk]; iCha++) {
              float Ylinfit = Par0 + Par1 * Xval[iDisk][iCha];
              hYdif[iDisk]->Fill(Ylinfit - Yval[iDisk][iCha]);
              Tres = Tval[iDisk][iCha] - Tcor[Ival[iDisk][iCha]] +
                     Yval[iDisk][iCha] / (cvel * costh) - t0;
              hTres[Ival[iDisk][iCha]]->Fill(Tres);

            } // Loop on readout channels
          } // Chi2 cut
        } // costh cut
      } // discr cut
    } // If enough good cells
  } // Loop on disks

  if (_diagLevel > 0)
    std::cout << "caloT0alig: end of event " << std::endl;

  return retval;

} // End of filter

// ===========================================================================
// End job:
// - Writing of temporary output file to be used from next iteration
// - Last iteration: Writing of calibration file to be uploaded in Condition DB
// ===========================================================================
void caloT0alig::endJob() {

  if (_diagLevel > 0)
    std::cout << "caloT0alig: Entering endJob" << std::endl;

  // Skip output writing when not in a minimization loop (i.e. filter only)
  if (_iteration == "filt") {
    if (_diagLevel > 0)
      std::cout << "Iteration " << _iteration << ": Skip output writing" << std::endl;
    return;
  }

  // Temporary output file for next iteration
  std::ofstream _fileTcor(_fileTcorName);
  if (!_fileTcor.is_open()) {
    throw cet::exception("caloT0alig")
        << "ERROR! Cannot open output file " << _fileTcorName << std::endl;
  }

  int Nevt;
  float Tmea, Tsig, Chi2, Tval;

  for (int iCha = 0; iCha < nROchan; iCha++) {
    TF1* gfit = new TF1("Gaussian", "gaus");
    hTres[iCha]->Fit(gfit, "q");
    Tmea = gfit->GetParameter(1);
    Tsig = gfit->GetParameter(2);
    Chi2 = gfit->GetChisquare();
    // TEMPORARY: MEAN/RMS instead of Gaussian fit params for MDC2020r CORSIKA production
    // Tmea = hTres[iCha]->GetMean();
    // Tsig = hTres[iCha]->GetRMS();
    Nevt = hTres[iCha]->GetEntries();
    Tval = Tmea + Tcor[iCha];

    _fileTcor << iCha << "   " << Tval << "   " << Tmea << "   " << Tsig << "   " << Chi2 << "   "
          << Nevt << std::endl;
  }
  _fileTcor.close();

  // LAST ITERATION:
  // Write output file for DB: rouId Tcor Tsigma Chi2 Nev
  // Link this calibration to the previous time calibration step from laser/cosmics
  // (unique DB identifier, table name + starting range, something else?)
  // What if different run numbers have different step-0 corrections?
  if (_iteration == "last") {
    if (_diagLevel > 0)
      std::cout << "Iteration " << _iteration << ": Writing output file for DB" << std::endl;
  }

} // End of endJob

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::caloT0alig);
