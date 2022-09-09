/*
 Author: Shihua Huang, Purdue University
 This program output pulse energy by analyzing pulse data read from a single txt file. The algorithm is defined in Mu2e-doc-30054. All the parameters are just temporary, needs to be optimized.
 The current data is taken from an oscilloscope, with built-in trigger. In the FPGA, we need to define our-own trigger. This part of code is acheived using a simple rising edge finder.
 The output should has pulse energy, pulse time and some quality factors(baseline shift, noise level,etc). The energy and time of the pulse is expected to be calculated separately in the FPGA.

Modified for Offline by A. Edmonds
*/

#include "Offline/STMReco/inc/PQAlg.hh"

namespace mu2e {

  int PQAlg::process_pulse(const mu2e::STMDigi& digi)
  {
    _energies.clear();
    _qualities.clear();
    _trigSamples.clear();
    //double sample_interval=2*digifac; //in nsec; not used
    int buffer_length=digi.adcs().size(); //224 points
    //    std::cout << "AE: buffer_length = " << buffer_length << std::endl;
    //comment the next line when not cutting the buffer;
    //buffer_length=cutfac;

    //!!test_y is to store the output energy!
    double test_y=0;

    //loop for a number of N_pulse pulse segments, which is the entire data file
    int i_0=0; //parameter used to increment through the first to the last data point(voltage value) in the data file; also marks the trigger value when a pulse rising edge is detected.
    int i_end=0; //marker to mark the 'end' position of one pulse segment
    //    int quality; //quality factor to check the pulse quality
    //    double x[buffer_length];
    double y[buffer_length];
    double peddle=0; // integrated pulse front pedestal value;
    double peddle2=0; // integrated pulse back pedestal value;
    double real_front=0;//the real front baseline value
    //    double slope;// the change in voltage w.r.t. time at a given time
    //    double i_trigger=0; // the pulse-start time by looking back from the trigger slope to the baseline, in one shnap-shot
    //    double peak=0; // the time of the CFD reference point, in one snap-shot
    //    double pulse_CFD_1=0; //CFD discrimination value1 of the pulse
    //    double pulse_CFD_2=0; //CFD discrimination value2 of the pulse

    //    double sample_interval=3.125;//4; //sample-interval is 4ns
    //*******************************
    //the main loop starts; once a trigger is found, it creates a 'pulse segement'(defined in Mu2e-doc-30054) and then 1) finds the 'real' pulse rise-up point (i_trigger) by walk back a little from the trigger point; 2) calculate the CFD based pulse time; 3) calculate pulse energy and 4) report pulse quality. In the current version result from 1) and 2) are not used.
    //*******************************
    j = digi.adcs().size();
    //    std::cout << "AE: " << i_0 << " < " << j << "-" << total_length << "-" << N_shift << " = " << j-total_length-N_shift <<"?" << std::endl;
    while(i_0<(j-total_length-N_shift)) //increment through the streaming data file's data points one by one; Make sure don't go exceed the end of data
      {
        // **************************************
        //pulse trigger detection and follow ups
        // **************************************
        if(//1024*1024*
           //           (digi.adcs().at(i_0)-digi.adcs().at(i_0+1))*(digi.adcs().at(i_0+1)-digi.adcs().at(i_0+2))<-2500)//0.25*sample_interval*sample_interval) //a trigger is set at i_0 whenever a rising edge is detected; a rising edge is detected if the multiply of the voltage difference between three consecutive data points,digi.adcs().at([i_0],digi.adcs().at(_[i_0+1] and digi.adcs().at([i_0+2] is larger than 1mV^2. (voltage expressed in mV)
           digi.adcs().at(i_0) < 3750)
          {
            //            std::cout << "AE: Triggered!" << std::endl;
            //            std::cout << "AE: (" << digi.adcs().at(i_0) << " - " << digi.adcs().at(i_0+1) << ")*(" << digi.adcs().at(i_0+1) << " - " << digi.adcs().at(i_0+2) << " = " << (digi.adcs().at(i_0)-digi.adcs().at(i_0+1))*(digi.adcs().at(i_0+1)-digi.adcs().at(i_0+2)) << " < -1000" << std::endl;//" << 0.25 << "*" << sample_interval << "^2 = " << 0.25*sample_interval*sample_interval << "?" << std::endl;
            //            std::cout << "AE: (" << digi.adcs().at(i_0) << ")" << std::endl;
            i_end=i_0+delta_Integration-delta_Walkback+Gap+delta_PB-1; //i_end is the end point of the pulse segment
            if (i_end < total_length) {
              i_0 = i_end;
              continue;
            }
            //            std::cout << "AE: i_end = " << i_end << std::endl;
            //parameter initilization
            peddle=0;peddle2=0;
            test_y=0;
            real_front=0;
            //            pulse_CFD_1=0;
            //            pulse_CFD_2=0;

            for (i=0;i<total_length;i++) //loop in each pulse segment, increment i
              {
                y[i]=digi.adcs().at(i+i_end-total_length+1);//*1000; //change pulse voltage in millie-volts units
                //                std::cout << "AE: y[i] = " << y[i] << std::endl;
              }

            //new section: 'real' rising-point-finder
            for(i=0;i<delta_PF;i++)
              {
                real_front=real_front+y[i]; //real front pedestal value
              }
            //trigger finding
            real_front=real_front/delta_PF;
            //            slope=1000*(digi.adcs().at(i_0)-digi.adcs().at(i_0+2))/2;//convert to millie-volts/bin
            //            i_trigger=delta_PF+delta_Walkback-(real_front-1000*digi.adcs().at(i_0+1))/slope;
            // slope=(y[i_trigger]-y[i_trigger+2])/2;
            // i_trigger=i_trigger-(real_front-y[i_trigger+1])/slope;


            //*************************************************
            //end of defining the trigger



            // **************************************
            //CFD pulse timing code
            // **************************************
            // for(i=0;i<total_length;i++) // loop in each snap-shot/pulse segment
            //   {
            //     pulse_CFD_1=pulse_CFD_2;
            //     pulse_CFD_2=
            //       a0*digi.adcs().at(i+i_end-total_length+1)-
            //       b0*digi.adcs().at(i+i_end-total_length+1+N_shift); //calculate the CFD_pulse value at position i;

            //     if((pulse_CFD_1+pulse_CFD_2)<0&&(pulse_CFD_1*pulse_CFD_2)<0)
            //       {peak=i-1+fabs(pulse_CFD_1)/fabs(pulse_CFD_2-pulse_CFD_1)-i_trigger;
            //         //fprintf(fp,"%f\n",peak); //CFD-pulse_start
            //         //fprintf(fp,"%f\n",peak);

            //       }//output pulse time at position i-1
            //     else if((pulse_CFD_1+pulse_CFD_2)>0&&(pulse_CFD_1*pulse_CFD_2)<0)
            //       {peak=i-1+fabs(pulse_CFD_1)/fabs(pulse_CFD_2-pulse_CFD_1)-i_trigger;
            //         //fprintf(fp,"%f\n",peak); //CFD-pulse_start
            //         // fprintf(fp,"%f\n",peak);

            //       } //output pulse time at position i
            //     else
            //       {; };

            //   }
            // **************************************
            //end of CFD pulse timing


            // **************************************
            //Simple pulse integration code
            // **************************************

            for(i=0;i<delta_PF;i++)
              {
                peddle=peddle+y[i]; //calculate the integrated front pedestal/baseline value
              }
            //            std::cout << "AE: peddle = " << peddle << std::endl;
            //            std::cout << "AE: peddle/delta_PF = " << peddle / delta_PF << std::endl;
            for(i=delta_PF;i<delta_PF+delta_Integration;i++)
              {
                test_y=test_y+(y[i]-3800);//(peddle/delta_PF));
                //pulse integration for energy(total charge) calculation
              }
            //            std::cout << "AE: test_y (already peddle subtracted) = " << test_y << std::endl;
                //            test_y=test_y-4*peddle; //subtract the baseline; the multiplicator 4 is because the integration length is 4 times longer than the front_pedestal
            //            std::cout << "AE: test_y (pre pC conversion) = " << test_y << std::endl;
            test_y *= -1.0/delta_Integration; // flip sign because negative pulse
            //            test_y=test_y*sample_interval/50;//convert to pC
            //***********************************
            //some diagnostic code; not used;
            // if(fabs(test_y)/128<24.5&&fabs(test_y)/128>24)
            //{//printf("%d\n",i_0);
            //
            //}

            //***********************************
            //some diagnostic code; not used;
            //if(/*i_0==43416*/i_0==85239)
            // { //printf("%d\n",i_0);
            //    for(i=0;i<total_length;i++)
            //  {//fprintf(fd,"%f\t%f\n",y[i],(y[i]-y[i+1])*(y[i+1]-y[i+2]));
            //    }
            // }
            // **************************************
            //end of simple pulse integration code
            // **************************************
            //pulse quality check
            // **************************************
            //            std::cout << "AE: test_y = " << test_y << std::endl;
            for(i=delta_PF+delta_Integration+Gap;i<total_length;i++)
              {
                peddle2=peddle2+y[i];//back pedestal value
              }
            if((peddle-peddle2<0)||(peddle-peddle2)*4*bound1>fabs(test_y)) //calculate the integrated front_back difference relatively large; May have pile-up
              {_qualities.push_back(2);
                _energies.push_back(test_y);
                _trigSamples.push_back(i_0);
                //fprintf(fp,"%f\t%d\t%d\n",test_y,i_0,quality);

              } //pile-up
            else if((peddle-peddle2)*4*bound2>fabs(test_y)) //front_back difference acceptable
              {_qualities.push_back(1);
                _energies.push_back(test_y);
                _trigSamples.push_back(i_0);
                //                fprintf(fp,"%f\t%d\t%d\n",test_y,i_0,quality); //okay quality, quality factor 1;
              }
            else //front_back difference is close to ideal
              {_qualities.push_back(0);
                _energies.push_back(test_y);
                _trigSamples.push_back(i_0);
                //                fprintf(fp,"%f\t%d\t%d\n",test_y,i_0,quality);
              } //relative good, qualify factor 0;


            // **************************************
            //end of pulse quality check


            i_0=i_end; //move i_0 to the end of pulse to perpare searching for the next trigger
          }

        // **************************************
        // end of the current trigger

        i_0=i_0+1; //increment i_0 and start searching for the next trigger
      }

    // **************************************
    // end of incrementing through the streaming data file
    //the main loop ends
    // **************************************

    //    fclose(fp);
    //fclose(fd);
    return 0;
  }
}
