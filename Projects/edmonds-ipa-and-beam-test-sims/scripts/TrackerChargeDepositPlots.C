#include "BaseRun.h"

#include "TH2.h"

TH2F* TrackerChargeDepositPlots(std::string id, std::vector<BaseRun*>& particle_runs, BaseRun* ce_run) {

  int n_years = BaseRun::n_years_running;
  int n_days = 365;
  int n_time_bins = n_years*52;

  double min_charge = 0;
  double max_charge = 10;
  double charge_width = 0.01;
  int n_charge_bins = (max_charge - min_charge) / charge_width;

  std::string histname = "hTotalChargeDeposit_" + id;
  TH2F* hTotalChargeDeposit = new TH2F(histname.c_str(), "", n_time_bins,0,n_time_bins, n_charge_bins,min_charge,max_charge);
  hTotalChargeDeposit->SetXTitle("Weeks of Running");
  hTotalChargeDeposit->SetYTitle("Straw Charge Deposit [C/cm]");
  hTotalChargeDeposit->SetZTitle("Number of Straw-cm");
  //  hTotalChargeDeposit->GetZaxis()->SetLabelOffset(0.00005);
  hTotalChargeDeposit->GetZaxis()->SetTitleOffset(0.75);
  hTotalChargeDeposit->SetStats(false);
  
  double charge_limit = 1; // C/cm
  histname = "hDeadStrawCms_" + id;
  TH1F* hDeadStrawCms = new TH1F(histname.c_str(), "", n_time_bins,0,n_time_bins);
  hDeadStrawCms->Sumw2();
  hDeadStrawCms->SetXTitle("Weeks of Running");
  hDeadStrawCms->SetYTitle("Fraction of Dead Straw-cm");

  histname = "hDeadStraws_" + id;
  TH1F* hDeadStraws = (TH1F*) hDeadStrawCms->Clone(histname.c_str());
  hDeadStraws->SetYTitle("Fraction of Dead Straws");

  TH3F* hConvEMinus_StrawCm_HitMap = ce_run->GetHitMap(); // get the conversion electron hitmape
  histname = "hFractionCEHitsLost_StrawCm_" + id;
  TH1F* hFractionCEHitsLost_StrawCm = (TH1F*) hDeadStrawCms->Clone(histname.c_str());
  hFractionCEHitsLost_StrawCm->SetYTitle("Fraction of CE Hits Lost (due to dead straw-cm)");

  TH2D* hConvEMinus_Straw_HitMap = (TH2D*) hConvEMinus_StrawCm_HitMap->Project3D("yz"); // get the conversion electron hitmap
  histname = "hFractionCEHitsLost_Straw_" + id;
  TH1F* hFractionCEHitsLost_Straw = (TH1F*) hDeadStraws->Clone(histname.c_str());
  hFractionCEHitsLost_Straw->SetYTitle("Fraction of CE Hits Lost (due to dead straws)");


  TAxis* straw_length_axis = (particle_runs.at(0))->GetHitMap()->GetXaxis(); // just so know what to loop over
      
  // Start looping through individual straws
  for (int i_device = 0; i_device < BaseRun::n_devices; ++i_device) {
    for (int i_sector = 0; i_sector < BaseRun::n_sectors; ++i_sector) {

      for (int i_straw = 0; i_straw < BaseRun::n_straws; ++i_straw) {
	int y_bin = 6*i_device + i_sector;
	int z_bin = i_straw;

	// Now loop along each cm of this straw
	bool straw_dead = false;
	int earliest_time = n_time_bins;
	for (int i_cm_bin = 1; i_cm_bin <= straw_length_axis->GetNbins(); ++i_cm_bin) {

	  // Loop through the different particle types
	  double charge_per_length = 0; // want to keep track of this for all particle types
	  for(std::vector<BaseRun*>::const_iterator i_particle_run = particle_runs.begin(); i_particle_run != particle_runs.end(); ++i_particle_run) {
	    TH3F* hHitMap = (*i_particle_run)->GetHitMap(); // this gives the number of hits per microbunch
	  
	    double hits_per_microbunch = hHitMap->GetBinContent(i_cm_bin, y_bin, z_bin);
	    charge_per_length += (hits_per_microbunch * (*i_particle_run)->GetChargeDepositPerHit() * BaseRun::n_microbunches_per_year * BaseRun::n_years_running); // over the whole experiment
	  }

	  // Now start looping through the time bins
	  for (int i_time_bin = 1; i_time_bin <= n_time_bins; ++i_time_bin) {
	    double time = hTotalChargeDeposit->GetXaxis()->GetBinLowEdge(i_time_bin);
	    
	    double charge_per_length_now = charge_per_length;
	    charge_per_length_now /= n_time_bins;
	    charge_per_length_now *= i_time_bin; // divide down to the time after this bin
	    
	    hTotalChargeDeposit->Fill(time, charge_per_length_now);
	    
	    if (charge_per_length_now >= charge_limit) {
	      hDeadStrawCms->Fill(time);
	      double fraction_ce_hits_lost_straw_cm = (hConvEMinus_StrawCm_HitMap->GetBinContent(i_cm_bin, y_bin, z_bin) * ce_run->GetNMicrobunchesSimulated()) / hConvEMinus_StrawCm_HitMap->GetEntries();
	      hFractionCEHitsLost_StrawCm->Fill(time, fraction_ce_hits_lost_straw_cm);

	      straw_dead = true;
	      if (i_time_bin < earliest_time) {
		earliest_time = i_time_bin;
	      }
	    }
	  }
	}

	// After we've gone across the whole straw, we know at what time it will die
	if (straw_dead) {
	  for (int i_time_bin = earliest_time; i_time_bin <= n_time_bins; ++i_time_bin) {
	    double time = hTotalChargeDeposit->GetXaxis()->GetBinLowEdge(i_time_bin);
	    hDeadStraws->Fill(time);
	    
	    double fraction_ce_hits_lost_straw = (hConvEMinus_Straw_HitMap->GetBinContent(y_bin, z_bin) * ce_run->GetNMicrobunchesSimulated()) / hConvEMinus_Straw_HitMap->GetEntries();
	    hFractionCEHitsLost_Straw->Fill(time, fraction_ce_hits_lost_straw);
	  }
	}
      }
    }
  }

  double total_number_of_straw_cm = 2.14031e6;
  hDeadStrawCms->Scale(1.0/total_number_of_straw_cm);
  hDeadStrawCms->SetLineWidth(2);
  hDeadStrawCms->SetStats(false);

  double total_number_of_straw = 23040;
  hDeadStraws->Scale(1.0/total_number_of_straw);
  hDeadStraws->SetLineWidth(2);
  hDeadStraws->SetStats(false);

  return hTotalChargeDeposit;
}
