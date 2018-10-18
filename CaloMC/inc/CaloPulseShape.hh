#ifndef CaloPulseShape_HH
#define CaloPulseShape_HH

#include <vector>
#include <string> 


namespace mu2e {


   class CaloPulseShape {

       public:

          CaloPulseShape(double digiSampling, int pulseIntegralSteps);
          ~CaloPulseShape() {};

          void buildShapes();
          
	  const std::vector< std::vector<double> >& pulseDigitized()      const  {return pulseDigitized_;}
          const std::vector<double>&                pulseDigitized(int i) const  {return pulseDigitized_.at(i);}
          const void                                printShape()          const;


      private:
      

	 double digiSampling_;
	 int    pulseIntegralSteps_;	 
	 std::vector< std::vector<double> > pulseDigitized_;
	 
   };

}
#endif
