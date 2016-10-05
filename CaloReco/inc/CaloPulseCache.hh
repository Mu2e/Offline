#ifndef CaloPulseCache_HH
#define CaloPulseCache_HH

#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include <vector>
#include <string> 


namespace mu2e {


   class CaloPulseCache {

       public:

          CaloPulseCache();
          ~CaloPulseCache() {};

	  void   initialize();
          double evaluate(double x);

          const std::vector<double>&   cache()      {return cache_;}
          double                       cache(int i) {return cache_.at(i);}
          double                       cacheSize()  {return cacheSize_;}
          double                       deltaT()     {return deltaT_;}
          double                       factor()     {return factor_;}
          double                       step()       {return step_;}


      private:

	 std::vector<double> cache_;
	 double              cacheSize_;
	 double              deltaT_;
	 double              factor_;
	 double              step_;
   };

}
#endif
