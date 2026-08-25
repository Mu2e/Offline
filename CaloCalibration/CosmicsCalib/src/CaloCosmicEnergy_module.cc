// =============================================================================
// extrapolate energy, position and time for offline calorimeter calibration
// 12-Dec-2023 - P.Fedeli & S.Giovannella
//
// Starting code for retrieving calo information:
// CaloEnergy & caloT0alig module from P.Fedeli & S.Giovannella
//
// =============================================================================

// includes

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/CrystalMapper.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// ROOT includes
#include "TF1.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"

namespace mu2e {

class CaloCosmicEnergy : public art::EDAnalyzer {

public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment(""), 0};
    fhicl::Atom<int> CutNCryHit{Name("CutNCryHit"), Comment(""), 3};
    fhicl::Atom<float> CutEnergyDep{Name("CutEnergyDep"), Comment(""), 10.};
    fhicl::Atom<float> CutChi2Norm{Name("CutChi2Norm"), Comment(""), 2.5};
    fhicl::Atom<art::InputTag> CaloClusterTag{
        Name("CaloClusterTag"), Comment("Tag for Calorimeter cluster collection"), art::InputTag()};
    fhicl::Atom<art::InputTag> CaloHitTag{
        Name("CaloHitTag"), Comment("Tag for Calorimeter hit collection"), art::InputTag()};
    fhicl::Atom<std::string> OutCalibFile{Name("OutCalibFile"),
                                          Comment("Name for output .dat MIP calibration file"),
                                          "calib_parameters.dat"};
    fhicl::Atom<std::string> OutCalibFileAsym{Name("OutCalibFileAsym"),
                                              Comment("Name for output .dat MIP asymmetry file"),
                                              "npe_noise.dat"};
  };

  explicit CaloCosmicEnergy(const art::EDAnalyzer::Table<Config>& config);
  virtual ~CaloCosmicEnergy() {};
  void beginJob() override;
  void beginRun(art::Run const& run) override;

  void analyze(art::Event const& event) override;

  void endJob() override;

  static double langaufun(Double_t* x, Double_t* par);

  float findpath(float m, float q, float x, float y);

private:
  art::ServiceHandle<art::TFileService> tfs;
  art::TFileDirectory tfdir = tfs->mkdir("All_tracks");
  art::TFileDirectory tfdirv = tfs->mkdir("Vertical_tracks");
  art::TFileDirectory tfdird = tfs->mkdir("Diagonal_tracks");
  art::TFileDirectory tfdirfp = tfs->mkdir("All_normalized_tracks");
  art::TFileDirectory tfdir_res = tfs->mkdir("MPV_distributions");
  art::TFileDirectory tfdir_cryalr = tfs->mkdir("Asymmetry_lr_Crystals");
  art::TFileDirectory tfdir_alr = tfs->mkdir("Analysis_asymmetry_lr");

  // calo parameters
  static constexpr int nSiPMs = CaloConst::_nSiPMPerCrystal;
  static constexpr int nCrystals = CaloConst::_nCrystalPerDisk;
  static constexpr int nDisks = CaloConst::_nDisk;
  static constexpr int nROchan = nDisks * nCrystals * nSiPMs;
  const Calorimeter* cal;          // For calorimeter geometry
  float cryDim;                    // Total crystal dimension
  float MaxDxVertical;             // Max Dx value for vertical tracks
  static constexpr int Erange = 5; // width of energy range in MeV for the sigma asymmetry variable
  static constexpr int Ebin = 11;  // number of ranges for the sigma asymmetry variable

  // fcl parameters
  int _diagLevel;
  int CutNCryHit;     //>
  float CutEnergyDep; // >
  float CutChi2Norm;  // <
  art::ProductToken<CaloClusterCollection> _caloClusterToken;
  art::ProductToken<CaloHitCollection> _caloHitToken;

  TH1F* hSiPM[nDisks * nCrystals][nSiPMs] = {nullptr};
  TH1F* hSiPMv[nDisks * nCrystals][nSiPMs] = {nullptr};
  TH1F* hSiPMd[nDisks * nCrystals][nSiPMs] = {nullptr};
  TH1F* hSiPMfp[nDisks * nCrystals][nSiPMs] = {nullptr};

  TH1F* MPV = nullptr;
  TH1F* MPVv = nullptr;
  TH1F* MPVd = nullptr;
  TH1F* MPVfp = nullptr;

  TH1F* LR_all = nullptr;
  TH2F* ALR_all = nullptr;
  TH1F* LR[Ebin] = {nullptr};
  TH1F* ALR[Ebin] = {nullptr};
  TH1F* CryALR[nDisks * nCrystals][Ebin] = {nullptr};

  TH1F* CryNpe = nullptr;
  TH1F* CrySigNoise = nullptr;
  TH1F* CryNpeChi2 = nullptr;

  float Energy_band[Ebin], counter_energy_band[Ebin], Cry_Energy_band[nDisks * nCrystals][Ebin],
      Cry_counter_energy_band[nDisks * nCrystals][Ebin];
  TF1* gaussianFit = new TF1("gaussianFit", "gaus");

  std::ofstream _outputfile;
  std::ofstream _outputfileAsym;

  // helper variables
  TF1* myPoly = new TF1("myPoly", "pol1", -600., 600.);    // polinomial to fit the track
  TF1* mylang = new TF1("mylang", langaufun, 15., 50., 4); // langaus function

}; // end class

CaloCosmicEnergy::CaloCosmicEnergy(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer(config), _diagLevel(config().diagLevel()), CutNCryHit(config().CutNCryHit()),
    CutEnergyDep(config().CutEnergyDep()), CutChi2Norm(config().CutChi2Norm()),
    _caloClusterToken(consumes<CaloClusterCollection>(config().CaloClusterTag())),
    _caloHitToken(consumes<CaloHitCollection>(config().CaloHitTag())),
    _outputfile(config().OutCalibFile()), _outputfileAsym(config().OutCalibFileAsym()) {

  if (!_outputfile.is_open()) {
    throw cet::exception("CaloCosmicEnergy")
        << "ERROR! Cannot open output file " << config().OutCalibFile() << std::endl;
  }
  if (!_outputfileAsym.is_open()) {
    throw cet::exception("CaloCosmicEnergy")
        << "ERROR! Cannot open output file " << config().OutCalibFileAsym() << std::endl;
  }
}

void CaloCosmicEnergy::beginRun(art::Run const& run) {

  art::ServiceHandle<GeometryService> geom;
  if (!geom->hasElement<Calorimeter>()) {
    throw cet::exception("NO-CALO-GEOM") << "CaloEnergy: Calorimeter geometry not found";
  }

  GeomHandle<Calorimeter> ch;
  cal = ch.get();
  cryDim = cal->G4Info().get<double>("crystalXYLength") +
           2. * cal->G4Info().get<double>("wrapperThickness");
  MaxDxVertical = cryDim * 1.1;
}

void CaloCosmicEnergy::beginJob() {
  if (_diagLevel > 0)
    std::cout << "CaloCosmicEnergy: Entering beginJob" << std::endl;

  for (int iCry = 0; iCry < nDisks * nCrystals; iCry++) {
    for (int iSiPM = 0; iSiPM < nSiPMs; iSiPM++) {
      hSiPM[iCry][iSiPM] = tfdir.make<TH1F>(
          Form("Crystal_%i_SiPM_%i", iCry, iSiPM),
          Form("Energy deposited in SiPM %i crystal %d", iSiPM, iCry), 70, 15., 50.);
      hSiPMv[iCry][iSiPM] = tfdirv.make<TH1F>(
          Form("Crystal_%i_SiPM_%i_v", iCry, iSiPM),
          Form("Energy deposited in SiPM %i crystal %d only vertical tracks", iSiPM, iCry), 70, 15.,
          50.);
      hSiPMd[iCry][iSiPM] = tfdird.make<TH1F>(
          Form("Crystal_%i_SiPM_%i_d", iCry, iSiPM),
          Form("Energy deposited in SiPM %i crystal %d only diagonal tracks", iSiPM, iCry), 70, 15.,
          50.);
      hSiPMfp[iCry][iSiPM] = tfdirfp.make<TH1F>(
          Form("Crystal_%i_SiPM_%i_fp", iCry, iSiPM),
          Form("Energy deposited in SiPM %i crystal %d normalized track length", iSiPM, iCry), 70,
          15., 50.);
    }
    for (int ibin = 0; ibin < Ebin; ibin++) {
      CryALR[iCry][ibin] = new TH1F(Form("Asymmetry_Crystal_%i_band_%i", iCry, ibin),
                                    Form("Asymmetry of crystal %i in band [%i, %i)", iCry,
                                         (ibin + 1) * Erange, (ibin + 1) * Erange + Erange),
                                    200, -1., 1.);
    }
  }
  MPV = tfdir_res.make<TH1F>("MPV_distr", "MPV distribution", 400, 10., 30.);
  MPVv =
      tfdir_res.make<TH1F>("MPV_v_distr", "MPV distribution only vertical tracks", 400, 10., 30.);
  MPVd =
      tfdir_res.make<TH1F>("MPV_d_distr", "MPV distribution only diagonal tracks", 400, 10., 30.);
  MPVfp =
      tfdir_res.make<TH1F>("MPV_fp_distr", "MPV distribution all normalized tracks", 400, 10., 30.);
  ALR_all =
      tfdir_alr.make<TH2F>("ALR", "Asymmetry left right all energies", 600, 0., 600., 100, -1., 1.);
  LR_all =
      tfdir_alr.make<TH1F>("ALR_All_Energies", "Asymmetry Left Right all energies", 200, -1, 1.);
  for (int i = 0; i < Ebin; i++) {
    LR[i] = tfdir_alr.make<TH1F>(
        Form("LR_ene_band_%i", i),
        Form("Ratio Left Right Energy band [%i - %i) MeV", i * Erange, (i + 1) * Erange), 100, 0.5,
        1.5);
    ALR[i] = tfdir_alr.make<TH1F>(
        Form("ALR_ene_band_%i", i),
        Form("Asymmetry Left-Right Energy band [%i - %i) MeV", i * Erange, (i + 1) * Erange), 200,
        -1., 1.);
  }
} // end begin job

void CaloCosmicEnergy::analyze(art::Event const& event) {

  // art::Handle<CaloClusterCollection> caloClustersHandle;
  // event.getByLabel(_caloClusterTag, caloClustersHandle);
  // const CaloClusterCollection& caloClusters(*caloClustersHandle);
  auto const& caloClusters = event.getProduct(_caloClusterToken);
  int nCluster = caloClusters.size();

  for (int iClu = 0; iClu < nCluster; iClu++) {
    int ncry = 0;     // number of crystals passing the selection
    float Dx = 1400.; // to calculate the spread alond x and y axis
    float Dy = 1400.;
    if (_diagLevel > 0) {
      std::cout << "iClu: " << iClu << std::endl;
    }

    // Creating arrays to store the information for crystals passing the selection
    std::array<float, nCrystals> PosX;
    std::array<float, nCrystals> PosY;
    std::array<float, nCrystals> ErrPos;
    ErrPos.fill(9.81f); // Error for a uniform distribution: 34/sqrt(12) ( 34 = crystal dimension )
    std::array<int, nCrystals> IDs;
    std::array<float, nCrystals> path;
    std::array<int, nCrystals> whichHit;

    /*
    float PosY    [nCrystals] = {0.0f};
    float ErrPos  [nCrystals] = {9.81f}; // Error for a uniform distribution: 34/sqrt(12) ( 34 =
    crystal dimension ) int   IDs     [nCrystals] = {0}; float path    [nCrystals] = {0.0f}; int
    whichHit[nCrystals] = {0};
    */

    // Loop over crystals
    for (long unsigned int iCry = 0; iCry < caloClusters[iClu].caloHitsPtrVector().size(); iCry++) {
      // cut on energy of crystal in cluster
      if (caloClusters[iClu].caloHitsPtrVector()[iCry].get()->energyDep() > CutEnergyDep) {
        // counter n crystal above 10 MeV
        PosX[ncry] =
            cal->mu2eToDiskFF(
                   caloClusters[iClu].diskID(),
                   cal->crystal(caloClusters[iClu].caloHitsPtrVector()[iCry].get()->crystalID())
                       .position())
                .getX();
        PosY[ncry] =
            cal->mu2eToDiskFF(
                   caloClusters[iClu].diskID(),
                   cal->crystal(caloClusters[iClu].caloHitsPtrVector()[iCry].get()->crystalID())
                       .position())
                .getY();
        IDs[ncry] = caloClusters[iClu].caloHitsPtrVector()[iCry].get()->crystalID();
        whichHit[ncry] = iCry;
        ncry++;
      }

    } // end loop iCry in cluster

    // cut on ncrystalhit
    if (ncry > CutNCryHit) {

      float min_x = PosX[0];
      float max_x = PosX[0]; // sbsitture xlictyposx->posX

      // bubble sort to find Xmax and Xmin
      for (int h = 1; h < ncry; h++) {
        if (PosX[h] > max_x) {
          max_x = PosX[h];
          // cout << "Max: " << max << endl;
        }
        if (PosX[h] < min_x) {
          min_x = PosX[h];
          // cout << "Min: " << min << endl;
        }
        Dx = abs(max_x - min_x);
      }

      float min_y = PosY[0];
      float max_y = PosY[0];

      // bubble sort to find Xmax and Xmin
      for (int h = 1; h < ncry; h++) {
        if (PosY[h] > max_y) {
          max_y = PosX[h];
          // cout << "Max: " << max << endl;
        }
        if (PosY[h] < min_y) {
          min_y = PosY[h];
          // cout << "Min: " << min << endl;
        }
        Dy = abs(max_y - min_y);
      }

      if (_diagLevel > 0) {
        std::cout << "dx: " << Dx << " dy: " << Dy << std::endl;
      }

      // geometry track
      TGraph* gfirst = new TGraph(ncry, PosX.data(), PosY.data());
      auto gfit = std::make_unique<TGraphErrors>(ncry, PosX.data(), PosY.data(), ErrPos.data(),
                                                 ErrPos.data());
      gfit->SetMarkerSize(0.9);
      gfit->SetMarkerStyle(8);

      float chi2norm = 0.;

      TFitResultPtr fitresult = gfirst->Fit("myPoly", "SQ");
      TFitResultPtr fit = gfit->Fit("myPoly", "SQ");
      // 2 fits help convergence

      chi2norm = fit->Chi2() / fit->Ndf();

      // begin findpath

      float TrkSlope = -1;
      TrkSlope = myPoly->GetParameter(1);
      if (Dx < 33) {     // if vertical tracks fit could be wrong
        TrkSlope = 144.; // ca 89.6 deg
      }

      for (int h = 0; h < ncry; h++) {
        if (_diagLevel > 0) {
          std::cout << "using m:" << TrkSlope << " q: " << myPoly->GetParameter(0)
                    << " x: " << PosX[h] << " y: " << PosY[h] << " ene: "
                    << caloClusters[iClu].caloHitsPtrVector()[whichHit[h]].get()->energyDep()
                    << std::endl;
        }
        path[h] = findpath(TrkSlope, myPoly->GetParameter(0), PosX[h], PosY[h]);
      } // end findpath

      // only verical tracks
      if (Dx < MaxDxVertical) {
        for (int iCry = 0; iCry < ncry; iCry++) {
          for (int iSiPM = 0; iSiPM < nSiPMs; iSiPM++) {
            hSiPMv[IDs[iCry]][iSiPM]->Fill(caloClusters[iClu]
                                               .caloHitsPtrVector()[whichHit[iCry]]
                                               .get()
                                               ->recoCaloDigis()[iSiPM]
                                               .get()
                                               ->energyDep());
          }
        }
      }

      // cut on diagonal tracks
      if (chi2norm < 0.01) {
        for (int iCry = 0; iCry < ncry; iCry++) {
          for (int iSiPM = 0; iSiPM < nSiPMs; iSiPM++) {
            hSiPMd[IDs[iCry]][iSiPM]->Fill(caloClusters[iClu]
                                               .caloHitsPtrVector()[whichHit[iCry]]
                                               .get()
                                               ->recoCaloDigis()[iSiPM]
                                               .get()
                                               ->energyDep());
          }
        }
      }

      // cut on vertical tracks + good fit chi2 -> all tracks
      if (Dx < 35 || chi2norm < CutChi2Norm) {
        for (int iCry = 0; iCry < ncry; iCry++) {
          for (int iSiPM = 0; iSiPM < nSiPMs; iSiPM++) {
            hSiPM[IDs[iCry]][iSiPM]->Fill(caloClusters[iClu]
                                              .caloHitsPtrVector()[whichHit[iCry]]
                                              .get()
                                              ->recoCaloDigis()[iSiPM]
                                              .get()
                                              ->energyDep());
            // normalized tracks by fidpaths.h
            if (path[iCry] > 0){
              hSiPMfp[IDs[iCry]][iSiPM]->Fill(caloClusters[iClu]
                                                .caloHitsPtrVector()[whichHit[iCry]]
                                                .get()
                                                ->recoCaloDigis()[iSiPM]
                                                .get()
                                                ->energyDep() *
                                            cryDim / path[iCry]);
            }
          }
        }
      } // chi2 cut

      // asymmetry left-right and ratio left/right->  we use only trigger selection
      for (long unsigned int iCry = 0; iCry < caloClusters[iClu].caloHitsPtrVector().size();
           iCry++) {
        float left_sipm = caloClusters[iClu]
                              .caloHitsPtrVector()[iCry]
                              .get()
                              ->recoCaloDigis()[0]
                              .get()
                              ->energyDep();
        float right_sipm = caloClusters[iClu]
                               .caloHitsPtrVector()[iCry]
                               .get()
                               ->recoCaloDigis()[1]
                               .get()
                               ->energyDep();
        float sipm_mean_e = caloClusters[iClu].caloHitsPtrVector()[iCry].get()->energyDep();
        int crystal_id = caloClusters[iClu].caloHitsPtrVector()[iCry].get()->crystalID();
        LR_all->Fill((left_sipm - right_sipm) / (left_sipm + right_sipm));
        ALR_all->Fill(sipm_mean_e, (left_sipm - right_sipm) / (left_sipm + right_sipm));
        if (sipm_mean_e >= 0. && sipm_mean_e < Ebin * Erange) {
          int whichband = sipm_mean_e / Erange;
          Energy_band[whichband] += sipm_mean_e;
          counter_energy_band[whichband]++;
          LR[whichband]->Fill(left_sipm / right_sipm);
          ALR[whichband]->Fill((left_sipm - right_sipm) / (left_sipm + right_sipm));
          CryALR[crystal_id][whichband]->Fill((left_sipm - right_sipm) / (left_sipm + right_sipm));
          Cry_Energy_band[crystal_id][whichband] += sipm_mean_e;
          Cry_counter_energy_band[crystal_id][whichband]++;
        }
      }

      delete gfirst;

    } // end ncrystal cut

  } // end loop iClu in ncluster
} // end enalyze

void CaloCosmicEnergy::endJob() {
  std::cout << "CaloCosmicEnergy: entering endjob" << std::endl;
  TGraphErrors *NpeLR, *NpeALR, *sigmaELR, *sigmaEALR, *CryNpeALR[nDisks * nCrystals];

  // fit params
  int redo[nDisks * nCrystals][nSiPMs];
  int ibin = 0;
  float Xmean = 0;
  float Xsigma = 0;
  float params[nDisks * nCrystals][nSiPMs][5];

  for (int k = 0; k < nDisks * nCrystals; k++) {
    for (int j = 0; j < nSiPMs; j++) {
      redo[k][j] = 0;
      for (int i = 0; i < 5; i++) {
        params[k][j][i] = -1;
      }
    }
  }
  for (int iCry = 0; iCry < nDisks * nCrystals; iCry++) {
    for (int iSiPM = 0; iSiPM < nSiPMs; iSiPM++) {

      // set initial parameters

      ibin = hSiPM[iCry][iSiPM]->GetMaximumBin();
      Xmean = hSiPM[iCry][iSiPM]->GetBinCenter(ibin);
      Xsigma = hSiPM[iCry][iSiPM]->GetRMS();

      mylang->SetParameters(Xsigma, Xmean, hSiPM[iCry][iSiPM]->GetMaximum(), 1);
      hSiPM[iCry][iSiPM]->Fit("mylang", "Q", "", 10, 50);

      if (_diagLevel > 0) {
        std::cout << "fit complete " << std::endl;
      }
      if ((mylang->GetChisquare() / mylang->GetNDF()) < 1.8 &&
          (hSiPM[iCry][iSiPM]->GetMean() - mylang->GetParameter(1)) < 10.) {
        if (_diagLevel > 0) {
          std::cout << "Crystal: " << iCry << " SiPM: " << iSiPM
                    << " MPV: " << mylang->GetParameter(1) << std::endl;
        }

        MPV->Fill(mylang->GetParameter(1));

        params[iCry][iSiPM][0] = mylang->GetParameter(1);
        params[iCry][iSiPM][1] = mylang->GetParError(1);
        params[iCry][iSiPM][2] = mylang->GetParameter(0);
        params[iCry][iSiPM][3] = mylang->GetParError(0);
        params[iCry][iSiPM][4] = mylang->GetChisquare() / mylang->GetNDF();

      } else {
        redo[iCry][iSiPM] = 1;
      }

    } // end sipms loop
  } // end crystal loop
  if (_diagLevel > 0) {
    std::cout << " I am redoing fits" << std::endl;
  }

  // redo wrong fit
  for (int iCry = 0; iCry < nDisks * nCrystals; iCry++) {
    for (int iSiPM = 0; iSiPM < nSiPMs; iSiPM++) {
      if (redo[iCry][iSiPM] == 1) {
        ibin = hSiPM[iCry][iSiPM]->GetMaximumBin();
        Xmean = hSiPM[iCry][iSiPM]->GetBinCenter(ibin);
        Xsigma = hSiPM[iCry][iSiPM]->GetRMS();

        // set wider parameter
        mylang->SetParameters(Xsigma * 1.2, Xmean, hSiPM[iCry][iSiPM]->GetMaximum(), 1.1);

        hSiPM[iCry][iSiPM]->Fit("mylang", "Q", "", 10, 50);

        MPV->Fill(mylang->GetParameter(1));

        params[iCry][iSiPM][0] = mylang->GetParameter(1);
        params[iCry][iSiPM][1] = mylang->GetParError(1);
        params[iCry][iSiPM][2] = mylang->GetParameter(0);
        params[iCry][iSiPM][3] = mylang->GetParError(0);
        params[iCry][iSiPM][4] = mylang->GetChisquare() / mylang->GetNDF();

        if (_diagLevel > 0) {
          std::cout << "Crystal: " << iCry << " SiPM: " << iSiPM
                    << " MPV: " << mylang->GetParameter(1) << std::endl;
        }
      }
    }
  }
  _outputfile << "#Cry #SiPM #MPV #MPV err #sigm #sigma err #chi2/ndf \n\n";
  for (int iCry = 0; iCry < nDisks * nCrystals; iCry++) {
    for (int iSiPM = 0; iSiPM < nSiPMs; iSiPM++) {
      _outputfile << Form("%i %i %.3f %.3f %.3f %.3f %.3f\n", iCry, iSiPM, params[iCry][iSiPM][0],
                          params[iCry][iSiPM][1], params[iCry][iSiPM][2], params[iCry][iSiPM][3],
                          params[iCry][iSiPM][4]);
    }
    _outputfile << "\n";
  }

  _outputfile.close();

  // fits for vertical diagonal and nromalized tracks
  for (int iCry = 0; iCry < nDisks * nCrystals; iCry++) {
    for (int iSiPM = 0; iSiPM < nSiPMs; iSiPM++) {

      // set initial parameters

      ibin = hSiPMv[iCry][iSiPM]->GetMaximumBin();
      Xmean = hSiPMv[iCry][iSiPM]->GetBinCenter(ibin);
      Xsigma = hSiPMv[iCry][iSiPM]->GetRMS();

      mylang->SetParameters(Xsigma, Xmean, hSiPMv[iCry][iSiPM]->GetMaximum(), 1);
      hSiPMv[iCry][iSiPM]->Fit("mylang", "Q", "", 10, 50);
      MPVv->Fill(mylang->GetParameter(1));

      ibin = hSiPMd[iCry][iSiPM]->GetMaximumBin();
      Xmean = hSiPMd[iCry][iSiPM]->GetBinCenter(ibin);
      Xsigma = hSiPMd[iCry][iSiPM]->GetRMS();

      mylang->SetParameters(Xsigma, Xmean, hSiPMd[iCry][iSiPM]->GetMaximum(), 1);
      hSiPMd[iCry][iSiPM]->Fit("mylang", "Q", "", 10, 50);
      MPVd->Fill(mylang->GetParameter(1));

      ibin = hSiPMfp[iCry][iSiPM]->GetMaximumBin();
      Xmean = hSiPMfp[iCry][iSiPM]->GetBinCenter(ibin);
      Xsigma = hSiPMfp[iCry][iSiPM]->GetRMS();

      mylang->SetParameters(Xsigma, Xmean, hSiPMfp[iCry][iSiPM]->GetMaximum(), 1);
      hSiPMfp[iCry][iSiPM]->Fit("mylang", "Q", "", 10, 50);
      MPVfp->Fill(mylang->GetParameter(1));
    }
  }
  MPV->Fit("gaus", "Q", "", 10., 30.);
  MPVv->Fit("gaus", "Q", "", 10., 30.);
  MPVd->Fit("gaus", "Q", "", 10., 30.);
  MPVfp->Fit("gaus", "Q", "", 10., 30.);

  // start Asymmetry L-R analysis

  LR_all->Fit("gaus", "Q", "", -1., 1.);

  // helper variables
  float sigmaLR[Ebin], sigmaLRerr[Ebin], sigmaALR[Ebin], sigmaALRerr[Ebin], NpeLRvar[Ebin],
      NpeALRvar[Ebin], eband_mean[Ebin], errEne[Ebin], errLR[Ebin], errALR[Ebin],
      CrySigmaALR[nDisks * nCrystals][Ebin], CrySigmaALRErr[nDisks * nCrystals][Ebin],
      Cryeband_mean[nDisks * nCrystals][Ebin], CryEneErr[nDisks * nCrystals][Ebin];

  // function to fit Asymmetry L-R(ALR) and L/R(LR)
  TF1* functional_form_ALR =
      new TF1("functional_form_ALR", "TMath::Sqrt(0.5*(1/([0]*x)+([1]/x)^2))", 0., 100.);
  TF1* functional_form_LR =
      new TF1("functional_form_LR", "TMath::Sqrt(2*(1/([0]*x))+([1]/x)^2)", 0., 100.);

  // peparing histogram for Npe e SigNoise distribution
  CryNpe = tfdir_alr.make<TH1F>("Npe_distr", "Distribution of Npe of all crystalss", 120, 0., 60.);
  CrySigNoise = tfdir_alr.make<TH1F>("SigNoise_distr", "Distribution of SigNoise of all crystals",
                                     80, 0., 0.4);
  CryNpeChi2 =
      tfdir_alr.make<TH1F>("CryNpeChi2", "Distribution of #chi2/ndof of all crystals", 100, 0., 50);

  // initializating string for npe-noise file dat
  _outputfileAsym << "#Cry #Npe/MeV #Npe/MeVerr #Noise #Noiseerr \n\n";

  // preparing variables to make tgraph for every crystal
  for (int icry = 0; icry < nDisks * nCrystals; icry++) {
    for (int ibin = 0; ibin < Ebin; ibin++) {
      CryALR[icry][ibin]->Fit("gaussianFit", "Q", "", -1., 1.);
      CrySigmaALR[icry][ibin] = gaussianFit->GetParameter(2);
      CrySigmaALRErr[icry][ibin] = gaussianFit->GetParError(2);

      // avoid infinities due low statistics
      if (Cry_counter_energy_band[icry][ibin] == 0) {
        Cry_counter_energy_band[icry][ibin] = 1;
      }
      Cryeband_mean[icry][ibin] = Cry_Energy_band[icry][ibin] / Cry_counter_energy_band[icry][ibin];
      CryEneErr[icry][ibin] = 1 / TMath::Sqrt(Cry_counter_energy_band[icry][ibin]);
    } // end ibin for

    // tgraph for single chanel
    CryNpeALR[icry] = tfdir_cryalr.make<TGraphErrors>(Ebin, Cryeband_mean[icry], CrySigmaALR[icry],
                                                      CryEneErr[icry], CrySigmaALRErr[icry]);
    CryNpeALR[icry]->SetMarkerStyle(8);
    CryNpeALR[icry]->SetMarkerSize(0.9);
    functional_form_ALR->SetParameters(30., 0.2);

    // fit single channel
    CryNpeALR[icry]->Fit("functional_form_ALR", "Q");

    // fill npe-nois file dat with crsystal; npe/mev +/- err; noise +/- err
    _outputfileAsym << Form("%i %.3f %.3f %.3f %.3f\n", icry, functional_form_ALR->GetParameter(0),
                            functional_form_ALR->GetParError(0),
                            functional_form_ALR->GetParameter(1),
                            functional_form_ALR->GetParError(1));
    // eventualy chi2/ndf cut
    if (functional_form_ALR->GetChisquare() / functional_form_ALR->GetNDF() > 0) {
      CryNpeALR[icry]->Write(Form("Cry_%i_NPe", icry));
      CryNpe->Fill(functional_form_ALR->GetParameter(0)); // fill Npe
      CryNpe->Fit("gaus", "Q");
      CrySigNoise->Fill(functional_form_ALR->GetParameter(1)); // fill SigNoise
      CrySigNoise->Fit("gaus", "Q");
    }
    CryNpeChi2->Fill(functional_form_ALR->GetChisquare() / functional_form_ALR->GetNDF());
  } // end icry for

  _outputfileAsym.close();
  for (int j = 0; j < Ebin; j++) {
    eband_mean[j] = Energy_band[j] / counter_energy_band[j];
    errEne[j] = 1 / TMath::Sqrt(counter_energy_band[j]);

    LR[j]->Fit("gaussianFit", "Q", "", 0.5, 1.5);
    sigmaLR[j] = gaussianFit->GetParameter(2);
    sigmaLRerr[j] = gaussianFit->GetParError(2);
    NpeLRvar[j] = 2 / gaussianFit->GetParameter(2) / gaussianFit->GetParameter(2) / eband_mean[j];
    errLR[j] = 4 * gaussianFit->GetParError(2) / gaussianFit->GetParameter(2) /
               gaussianFit->GetParameter(2) / gaussianFit->GetParameter(2) / eband_mean[j];

    ALR[j]->Fit("gaussianFit", "Q", "", -1., 1.);
    sigmaALR[j] = gaussianFit->GetParameter(2);
    sigmaALRerr[j] = gaussianFit->GetParError(2);
    NpeALRvar[j] =
        1 / gaussianFit->GetParameter(2) / gaussianFit->GetParameter(2) / eband_mean[j] / 2;
    errALR[j] = 0.5 * gaussianFit->GetParError(2) / gaussianFit->GetParameter(2) /
                gaussianFit->GetParameter(2) / gaussianFit->GetParameter(2) / eband_mean[j];
  }

  // analsis for all cristals togheter
  NpeLR = tfdir_alr.make<TGraphErrors>(Ebin, eband_mean, NpeLRvar, errEne, errLR);
  NpeALR = tfdir_alr.make<TGraphErrors>(Ebin, eband_mean, NpeALRvar, errEne, errALR);
  sigmaELR = tfdir_alr.make<TGraphErrors>(Ebin, eband_mean, sigmaLR, errEne, sigmaLRerr);
  sigmaEALR = tfdir_alr.make<TGraphErrors>(Ebin, eband_mean, sigmaALR, errEne, sigmaALRerr);

  NpeLR->SetMarkerSize(0.9);
  NpeLR->SetMarkerStyle(8);
  NpeALR->SetMarkerSize(0.9);
  NpeALR->SetMarkerStyle(8);
  sigmaELR->SetMarkerSize(0.9);
  sigmaELR->SetMarkerStyle(8);
  sigmaEALR->SetMarkerSize(0.9);
  sigmaEALR->SetMarkerStyle(8);
  NpeLR->Fit("pol0", "Q");
  NpeALR->Fit("pol0", "Q");
  functional_form_LR->SetParameters(30., 0.2);
  functional_form_ALR->SetParameters(30., 0.2);
  sigmaELR->Fit("functional_form_LR", "Q");
  sigmaEALR->Fit("functional_form_ALR", "Q");
  NpeLR->Write("NpeLR");
  NpeALR->Write("NpeALR");
  sigmaELR->Write("sigmaELR");
  sigmaEALR->Write("sigmaEALR");

  if (_diagLevel > 0) {
    std::cout << "[0]: " << functional_form_ALR->GetParameter(0)
              << "[1]: " << functional_form_ALR->GetParameter(1);
    std::cout << "[0]: " << functional_form_LR->GetParameter(0)
              << "[1]: " << functional_form_LR->GetParameter(1);
  }

  delete gaussianFit;
  delete functional_form_ALR;
  delete functional_form_LR;
  for (int icry = 0; icry < nDisks * nCrystals; icry++) {
    for (int ibin = 0; ibin < Ebin; ibin++) {
      delete CryALR[icry][ibin];
    }
  }

} // end endjob

double CaloCosmicEnergy::langaufun(double* x, double* par) {

  // Fit parameters:
  // par[0]=Width (scale) parameter of Landau density
  // par[1]=Most Probable (MP, location) parameter of Landau density
  // par[2]=Total area (integral -inf to inf, normalization constant)
  // par[3]=Width (sigma) of convoluted Gaussian function
  //
  // In the Landau distribution (represented by the CERNLIB approximation),
  // the maximum is located at x=-0.22278298 with the location parameter=0.
  // This shift is corrected within this function, so that the actual
  // maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
  Double_t mpshift = -0.22278298;      // Landau maximum location

  // Control constants
  Double_t np = 100.0; // number of convolution steps
  Double_t sc = 5.0;   // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow, xupp;
  Double_t step;
  Double_t i;

  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp - xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for (i = 1.0; i <= np / 2; i++) {
    xx = xlow + (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);

    xx = xupp - (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
} // end langaus

float CaloCosmicEnergy::findpath(float m, float q, float x, float y) {
  float xup = 0;
  float yright = 0;
  float path = -1;
  float xlow = 0;
  float yleft = 0;

  int diag = 0;
  float halfcrysize = cryDim / 2.;

  if (m != 0) {
    xup = ((y + halfcrysize) - q) / m;
    xlow = ((y - halfcrysize) - q) / m;
    yleft = m * (x - halfcrysize) + q;
    yright = m * (x + halfcrysize) + q;
  }
  if (xup < (x + halfcrysize) && xup > (x - halfcrysize)) { // hit the top of crystal
    if (xlow < (x + halfcrysize) && xlow > (x - halfcrysize)) {
      path = halfcrysize * 2 * TMath::Sqrt(1 / (m * m) + 1);
      if (diag == 1) {
        std::cout << "x0,y0: (" << xup << ";" << y + halfcrysize << ") x1,y1: (" << xlow << ";"
                  << y - halfcrysize << ")" << std::endl;
      }
    }
    if (xlow < (x - halfcrysize)) {
      path = TMath::Sqrt((xup - (x - halfcrysize)) * (xup - (x - halfcrysize)) +
                         ((y + halfcrysize) - yleft) * ((y + halfcrysize) - yleft));
      if (diag == 1) {
        std::cout << "Track exit from left" << std::endl;
        std::cout << "x0,y0: (" << xup << ";" << y + halfcrysize << ") x1,y1: (" << x - halfcrysize
                  << ";" << yleft << ")" << std::endl;
      }
    }
    if (xlow > (x + halfcrysize)) {
      path = TMath::Sqrt((xup - (x + halfcrysize)) * (xup - (x + halfcrysize)) +
                         ((y + halfcrysize) - yright) * ((y + halfcrysize) - yright));
      if (diag == 1) {
        std::cout << "x0,y0: (" << xup << ";" << y + halfcrysize << ") x1,y1: (" << x + halfcrysize
                  << ";" << yright << ")" << std::endl;
        std::cout << "Track exit from right" << std::endl;
      }
    }
  }
  // enter from left
  else if (yleft < (y + halfcrysize) && yleft > (y - halfcrysize)) {
    if (diag == 1) {
      std::cout << "track enter from left and ";
    }
    if (xlow < (x + halfcrysize) && xlow > (x - halfcrysize)) {
      if (diag == 1) {
        std::cout << " it exit bot " << std::endl;
        std::cout << "x0,y0: (" << x - halfcrysize << ";" << yleft << ") x1,y1: (" << xlow << ";"
                  << y - halfcrysize << ")" << std::endl;
      }
      path = TMath::Sqrt((x - halfcrysize - xlow) * (x - halfcrysize - xlow) +
                         (yleft - (y - halfcrysize)) * (yleft - (y - halfcrysize)));
    } else { // track exit from right
      if (diag == 1) {
        std::cout << " it exit from right" << std::endl;
        std::cout << "x0,y0: (" << x - halfcrysize << ";" << yleft << ") x1,y1: ("
                  << x + halfcrysize << ";" << yright << ")" << std::endl;
      }
      path = halfcrysize * 2 * TMath::Sqrt(1 + m * m);
    }
  }

  else if (yright < (y + halfcrysize) && yright > (y - halfcrysize)) {
    if (diag == 1) {
      std::cout << "track enter from right and ";
    }
    if (xlow < (x + halfcrysize) && xlow > (x - halfcrysize)) {
      if (diag == 1) {
        std::cout << " it exit bot " << std::endl;
        std::cout << "x0,y0: (" << x + halfcrysize << ";" << yright << ") x1,y1: (" << xlow << ";"
                  << y - halfcrysize << ")" << std::endl;
      }
      path = TMath::Sqrt((x + halfcrysize - xlow) * (x + halfcrysize - xlow) +
                         (yright - (y - halfcrysize)) * (yright - (y - halfcrysize)));
    }
  }

  else {
    if (diag == 1) {
      std::cout << "can't find a path" << std::endl;
    }
  }
  if (diag == 1) {
    std::cout << "path is: " << path << std::endl;
  }

  return path;
} // end findpath

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloCosmicEnergy)
