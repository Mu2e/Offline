#include "BaseRun.h"

#include "TH2.h"

TH2F* StrawChargeDepositOverTime(std::string id, BaseRun* particle_run) {

  int n_years = BaseRun::n_years_running;
  int n_days = 365;
  int n_time_bins = n_years*52;

  double min_charge = 0;
  double max_charge = 10;
  double charge_width = 0.01;
  int n_charge_bins = (max_charge - min_charge) / charge_width;

  std::string histname = "hStrawChargeDeposit_" + id;
  TH2F* hStrawChargeDeposit = new TH2F(histname.c_str(), "", n_time_bins,0,n_time_bins, n_charge_bins,min_charge,max_charge);
  hStrawChargeDeposit->SetXTitle("Weeks of Running");
  hStrawChargeDeposit->SetYTitle("Straw Charge Deposit [C/cm]");
  hStrawChargeDeposit->SetZTitle("Number of Straw-cm");
  //  hStrawChargeDeposit->GetZaxis()->SetLabelOffset(0.00005);
  hStrawChargeDeposit->GetZaxis()->SetTitleOffset(0.75);
  hStrawChargeDeposit->SetStats(false);
  

  TAxis* straw_length_axis = particle_run->GetHitMap()->GetXaxis(); // just so know what to loop over
      
  // Start looping through individual straws
  for (int i_device = 0; i_device < BaseRun::n_devices; ++i_device) {
    for (int i_sector = 0; i_sector < BaseRun::n_sectors; ++i_sector) {
      for (int i_straw = 0; i_straw < BaseRun::n_straws; ++i_straw) {
	int y_bin = 6*i_device + i_sector;
	int z_bin = i_straw;

	// Now loop along each cm of this straw
	for (int i_cm_bin = 1; i_cm_bin <= straw_length_axis->GetNbins(); ++i_cm_bin) {
	  TH3F* hHitMap = particle_run->GetHitMap(); // this gives the number of hits per microbunch
	  
	  double hits_per_microbunch = hHitMap->GetBinContent(i_cm_bin, y_bin, z_bin);
	  double charge_per_length = (hits_per_microbunch * particle_run->GetChargeDepositPerHit() * BaseRun::n_microbunches_per_year * BaseRun::n_years_running); // over the whole experiment

	  // Now start looping through the time bins
	  for (int i_time_bin = 1; i_time_bin <= n_time_bins; ++i_time_bin) {
	    double time = hStrawChargeDeposit->GetXaxis()->GetBinLowEdge(i_time_bin);

	    double charge_per_length_now = charge_per_length;
	    charge_per_length_now /= n_time_bins;
	    charge_per_length_now *= i_time_bin; // divide down to the time after this bin

	    hStrawChargeDeposit->Fill(time, charge_per_length_now);
	  }
	}
      }
    }
  }

  return hStrawChargeDeposit;
}
