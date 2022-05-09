///////////////////////////////////////
//
// Zero-suppression algorithm for STM
// Original authors: Claudia Alvarez-Garcia, Alex Kesharvarzi, and Mark Lancaster
// Adapated for Offline: Andy Edmonds
//
#include "Offline/STMReco/inc/ZPAlg.hh"

namespace mu2e {

  int16_t* ZPAlg::ZeroSuppress(const int16_t* ADC, unsigned long int n, std::vector<size_t>& starts, std::vector<size_t>& ends){
    //n number of elements to apply zero sup

    //----------------------- PEAKS SUPPRESSION ALGORITHM-----------------------------------------

    //Calculate the gradient for ADC values separated a distance window (in time would be window * tadc)
    int window=100;
    int16_t* gradient = new int16_t[n-window];
    //  double* time = new double[n-window];
    int16_t* suppressed_data = new int16_t[n-window];
    //Sampling time of ADC (microsec)
    const double tadc=1/(fADC*1e-6);

    //Calculate the gradient
    for(long unsigned int i=0;i<(n-window);i++){

      //time in microsec
      //time[i]=i*tadc;
      gradient[i]=ADC[i+window]-ADC[i];
    }


    //Initial values
    bool peak=false;
    int counter=0;
    int trigger=0;
    int triggercounter=0;


    //store 2 microseconds of data to the left of the trigger
    int prenumADCstored=int(tbefore/tadc);

    //store 20 microseconds of data to the right of the trigger
    int postnumADCstored=int(tafter/tadc);
    //---------------------------


    for(unsigned long int i=0;i<(n-window);i++){

      if(gradient[i]>threshold){peak=false;
        continue;}
      if((gradient[i]<threshold)&&(peak==true)){continue;}

      int triggerold=trigger;
      if((gradient[i]<threshold)&&(peak==false)){
        trigger=i;
        //if triggers found are closer than the window of stored data in time, rejected trigger, continue
        //it is necesary to convert the timing window to the index window: time[i]=i*tadc
        if(trigger-triggerold<int((tafter+tbefore)/tadc)){continue;}

        cout<<"Trigger number: "<<triggercounter<<": "<<trigger<<" Triggertime: "<<trigger*tadc<<endl;
        triggercounter++;

        for(int k=prenumADCstored; k>0;k--){suppressed_data[counter]=ADC[trigger-k];
          peak=true;
          counter++;}

        for(int j=0; j<postnumADCstored;j++){suppressed_data[counter]=ADC[trigger+j];
          peak=true;
          counter++;}

        starts.push_back(i-prenumADCstored);
        ends.push_back(i+postnumADCstored);
      }
    }

    std::cout <<"Number of triggers/peaks found: "<<triggercounter<<std::endl;
    std::cout<<"Number of elements in suppressed file: "<<counter<<std::endl;

    //Return array of suppressed data
    return suppressed_data;
  }
}
