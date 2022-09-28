#include "Offline/STMReco/inc/MWDAlg.hh"
#include <chrono>
#include <thread>

namespace mu2e {
  /*
  MWDAlg::MWDAlg(){
    M = 1000;
    L = 500;
    tau = 50000;  // in ns
    nsigma_cut = 9;
    thresholdgrad = -0.3;
    fADC = 320.0520833313; // in MHz
    cut_mode = 1;
    fixed_cut_parameter = -1000.0;
  }
  */
  MWDAlg::MWDAlg(double _M, double _L, double _tau, double _nsigma_cut, double _thresholdgrad, double _fADC, int _cut_mode, double _fixed_cut_parameter){
    M = _M;
    L = _L;
    tau = _tau;
    nsigma_cut = _nsigma_cut;
    thresholdgrad = _thresholdgrad;
    fADC = _fADC;
    cut_mode = _cut_mode;
    fixed_cut_parameter = _fixed_cut_parameter;
  }

  std::string MWDAlg::print() {
    std::stringstream ss;
    ss << "M             = " << M << "\n"
       << "L             = " << L << "\n"
       << "tau           = " << tau << "\n"
       << "nsigma_cut    = " << nsigma_cut << "\n"
       << "fixed_cut     = " << fixed_cut_parameter << "\n"
       << "cut_mode      = " << cut_mode << "\n"
       << "nsigma_cut    = " << nsigma_cut << "\n"
       << "thresholdgrad = " << thresholdgrad << "\n"
       << "fADC [MHz]    = " << fADC << "\n" ;
    return ss.str();
  }

  peaks* MWDAlg::find_peaks(double baseline_mean, double baseline_rms, double time_offset) {

    peaks* peak_data = new peaks();
    peak_data->npeaks = 0;

    if (cut_mode == 1)
      threshold_cut = baseline_mean - baseline_rms*nsigma_cut;
    else if (cut_mode == 2)
      threshold_cut = fixed_cut_parameter;
    else {
      std::cout << "ERROR:: find_peaks() : unknown cut mode .... " << cut_mode << std::endl;
      exit(-1);
    }

    //    if (cut_mode == 1)
      //      std::cout << "MWDAlg::find_peaks()---SIGMA CUT : baseline mean "<< baseline_mean << " rms = " << baseline_rms << " num_sigmas = " << nsigma_cut
      //                << " threshold cut value = " << threshold_cut << " and M = " << M << " L = " << L << " and MWDAlg cut mode = " << cut_mode << std::endl;
    //    if (cut_mode == 2)
      //      std::cout << "MWDAlg::find_peaks()---FIXED CUT : baseline mean "<< baseline_mean << " rms = " << baseline_rms
      //                << " threshold cut value = " << threshold_cut << " and M = " << M << " L = " << L << " and MWDAlg cut mode = " << cut_mode << std::endl;

    int n = nadc;
    double* e2 = new double[n];
    double* energy = new double[n];
    double* timepeak = new double[n];
    double* adc_time = new double[n];
    double auxlow = 0.0 ;
    int counterpeak = 0 ;
    double timeaux = 0;

    for(int i = 0; i < n; i++){
      adc_time[i] = time_offset + ((double) i)/( fADC  ) ; // in us
    }

    for( int i = M; i < n; i++){

      if (l[i] < threshold_cut){

        if ((l[i] < l[i-1]) && l[i] < auxlow){

          auxlow = l[i];
          timeaux = adc_time[i];
          e2[counterpeak] = auxlow;

        }
        else {
          continue;
        }
      }

      if (auxlow == 0) {
        continue;
      }
      else if (l[i] > threshold_cut){

        peak_data->npeaks++;
        energy[counterpeak] = e2[counterpeak] - baseline_mean;
        peak_data->peak_heights.push_back(energy[counterpeak]);
        peak_data->peak_times.push_back(timeaux); //us

        //std::cout << "Peak[" << counterpeak << "], Time: " << timeaux <<" us" << ", Energy (ADC counts): " << energy[counterpeak] << std::endl;
        auxlow=0.;
        counterpeak++;
      }
    }
    delete e2; delete energy; delete timepeak; delete adc_time;
    return peak_data;
  }

  std::vector<double> MWDAlg::calculate_baseline(){

    int k = M;
    int n = nadc;

    //std::cout << "l values in calculate_baseline.... " << l[0] << " " << l[1] << " " << l[2] << std::endl;

    double* gradient = new double[n];
    double* lvalues = new double[n];

    int ilv = 0;

    //Remove peaks and calculate MWDAlg baseline mean
    while (k < n){
      gradient[k] = l[k+1] - l[k];
      //std::cout << "AE: grad[" << k << "] = " << gradient[k] << ", l[" <<  k << "] = " << l[k] << std::endl;
      if(gradient[k] < thresholdgrad){
        k = k + (M+2*L);
        continue;
      }
      else {
        lvalues[ilv] = l[k];
        ilv++;
        k++;
      }
    }
    double mean = stats::mean(lvalues,ilv);
    double rms = 0;
    if (ilv > 1) {
      rms = stats::rms(lvalues,mean,ilv);
    }
    std::vector<double> result;
    // std::cout << " Baseline Mean " << mean << " RMS = " << rms << std::endl;
    result.push_back(mean);
    result.push_back(rms);

    delete lvalues;
    delete gradient;

    return result;
  }

  void MWDAlg::mwd_algorithm(data* adc_values){  // This fill the double[l] array ie sets the private l*
    int n = adc_values->nadc;
    nadc = n;
    //std::cout << "adc_values .... " << adc_values->adc[0] << " " << adc_values->adc[1] << " " << adc_values->adc[2] << std::endl;
    //std::cout << "n = " << n << std::endl;
    const double T0 = (1000.0/fADC); // in ns

    ////////////////////////////////     MWD Algorithm   //////////////////////////////////////////////////
    /*
    std::ofstream out1, out2, out3, out4;
    out1.open("MWDStepDecon.log", std::ios::out | std::ios::trunc);
    out1 << "Index_i," << "a_i" << std::endl;
    out2.open("MWDStepDiff.log", std::ios::out | std::ios::trunc);
    out2 << "Index_i," << "D_i" << std::endl;
    out3.open("MWDStepAvg.log", std::ios::out | std::ios::trunc);
    out3 << "Index_i," << "l_i" << std::endl;
    out4.open("MWDStepWaveform.log", std::ios::out | std::ios::trunc);
    out4 << "Index_i," << "adc_i" << std::endl;
    for (int i = 1; i < n; ++i)
      {
        out4 << i << ", " << adc_values->adc[i] << std::endl;
      }
    */
    //Deconvolution
    double* a = new double[n];
    double pedestal_sum = 0;
    for (int i = 0; i < 100; ++i)
      {
        pedestal_sum += adc_values->adc[i];
      }
    double pedestal = pedestal_sum/100;
    a[0] = adc_values->adc[0];
    for(int i=1; i<n; i++){
      a[i] = (adc_values->adc[i]-pedestal)-(1-(T0/tau))*(adc_values->adc[i-1]-pedestal) + a[i-1];
      //out1 << i << ", " << a[i] << std::endl;
    }
    //std::cout << "adc_values .... " << adc_values->adc[0] << " " << adc_values->adc[1] << " " << adc_values->adc[2] << std::endl;
    //std::cout << "a values .... " << a[0] << " " << a[1] << " " << a[2] << std::endl;

    //Differentiation
    double* D = new double[n];
    memcpy( D, a, M*sizeof(a) );

    for (int i = M; i < n; ++i) {
      D[i] = a[i] - a[i-M];
    }
    /*
    for(int i = 0; i < n; ++i)
      {
         out2 << i << ", " << D[i] << std::endl;
      }
    */
    //std::cout << "adc_values .... " << adc_values->adc[0] << " " << adc_values->adc[1] << " " << adc_values->adc[2] << std::endl;
    //std::cout << "D values[0..2] .... " << D[0] << " " << D[1] << " " << D[2] << std::endl;
    //std::cout << "D values[M..M+2].... " << D[M] << " " << D[M+1] << " " << D[M+2] << std::endl;
    delete a;
    //std::cout << "adc_values .... " << adc_values->adc[0] << " " << adc_values->adc[1] << " " << adc_values->adc[2] << std::endl;

    //Averaging
    l = new double[n];
    //std::cout << "l values before memcpy.... " << l[0] << " " << l[1] << " " << l[2] << std::endl;
    double sum = 0.;

    //std::cout << "size of l array = " << n << " nbyes being copied = " << (L-1)*sizeof(D) << std::endl;
    memcpy( l, D, (L-1)*sizeof(D) );
    //std::cout << "l values after memcpy.... " << l[0] << " " << l[1] << " " << l[2] << std::endl;
    for (int i = 0; i < L-1; ++i) {
      sum += D[i];
    }
    //std::cout << " sum = " << sum << std::endl;

    sum += D[L-1];
    //std::cout << " sum = " << sum << std::endl;
    l[L-1] = sum/L;
    //std::cout << " l[L-1] = " << l[L-1] << std::endl;

    for (int i = L; i < n; ++i) {
      sum += D[i]-D[i-L];
      l[i] = sum/L;
    }
    /*
    for (int i = 0; i < n; ++i)
      {
        out3 << i << ", " << l[i] << std::endl;
      }
    */
    //std::cout << " sum = " << sum << std::endl;
    //for (int i = 0; i < 10; i++) {
    //  std::cout << "l[" << i << "] = " << l[i] << std::endl;
    //}

    delete D;
  }
};
