//
// Contains data for the calorimeter class. A lot of information should move to a DB
//
// Original author B. Echenard
//
#ifndef CalorimeterGeom_CaloInfo_hh
#define CalorimeterGeom_CaloInfo_hh

#include <vector>

namespace mu2e {

    class CaloInfo {

       public:

           CaloInfo(): 
	      envelopeRadiusIn_(0),envelopeRadiusOut_(0),envelopeZ0_(0),envelopeZ1_(0),             
              crystalShift_(0),crystalXYLength_(0),crystalZLength_(0),crystalVolume_(0),
              wrapperThickness_(0),refractiveIndex_(0),crystalDecayTime_(0),
	      nROPerCrystal_(0),roXLength_(0),roYLength_(0),roZLength_(0),
              FEEXLength_(0),FEEYLength_(0),FEEZLength_(0),FEEBoxThickness_(0),
              BPHoleXLength_(0),BPHoleYLength_(0),BPHoleZLenght_(0),stripThickness_(0),stripLengthY_(0),coolBPPipeRadius_(0),
              outerRingEdgeZLength_(0),outerRingEdgeRLength_(0),caseThicknessIn_(0),caseThicknessOut_(0),
              FPStyrofoamZLength_(0),FPCarbonZLength_(0),coolFPPipeRadius_(0),
	      nPipes_(0),pipeRadius_(0), pipeThickness_(0),pipeTorRadius_(),pipeSeparation_(0),              
              nCrate_(0),nBoard_(0),nCrateBeforeSpace_(0),crateXLength_(0),crateYLength_(0),crateZLength_(0),
              crateFShieldThick_(0),crateBShieldThick_(0),crateTThick_(0),crateSThick_(0),crateFShieldYLength_(0),
              crateFShieldDeltaZ_(0),crateRadiusIn_(0),cratephi0_(0),crateDeltaPhi_(0),
              radiatorThickness_(0),activeStripThickness_(0),passiveStripThickness_(0)
           {}
	     
           ~CaloInfo() {}
           
	   int crystalByRO(int roid)          const  {return (roid/nROPerCrystal_);}
	   int ROBaseByCrystal(int crystalId) const  {return (crystalId*nROPerCrystal_);}


           void envelopeRadiusIn(double value)     {envelopeRadiusIn_ = value;}
           void envelopeRadiusOut(double value)    {envelopeRadiusOut_ = value;}
           void envelopeZ0(double value)           {envelopeZ0_ = value;}
           void envelopeZ1(double value)           {envelopeZ1_ = value;} 
           
           void crystalShift(bool value)           {crystalShift_ = value;}
           void crystalZLength(double value)       {crystalZLength_ = value;}
           void crystalXYLength(double value)      {crystalXYLength_ = value;}
           void crystalVolume(double value)        {crystalVolume_ = value;}
	   void wrapperThickness(double value)     {wrapperThickness_ = value;}
	   void refractiveIndex(double value)      {refractiveIndex_ = value;}
	   void crystalDecayTime(double value)     {crystalDecayTime_ = value;}
           
	   void nROPerCrystal(int value)           {nROPerCrystal_ = value;}
           void roXLength(double value)            {roXLength_ = value;}
           void roYLength(double value)            {roYLength_ = value;}
           void roZLength(double value)            {roZLength_ = value;}
           void FEEXLength(double value)           {FEEXLength_ = value;}
           void FEEYLength(double value)           {FEEYLength_ = value;}
           void FEEZLength(double value)           {FEEZLength_ = value;}
           void FEEBoxThickness(double value)      {FEEBoxThickness_=value;}
           void BPHoleXLength(double value)        {BPHoleXLength_ = value;}
           void BPHoleYLength(double value)        {BPHoleYLength_ = value;}
           void BPHoleZLength(double value)        {BPHoleZLenght_ = value;}
           void stripThickness(double value)       {stripThickness_=value;}
           void stripYLength(double value)         {stripLengthY_ = value;}
           void coolBPPipeRadius(double value)     {coolBPPipeRadius_ = value;}
           
           void outerRingEdgeZLength(double value) {outerRingEdgeZLength_ = value;}
           void outerRingEdgeRLength(double value) {outerRingEdgeRLength_ = value;}
           void caseThicknessIn(double value)      {caseThicknessIn_ = value;}
           void caseThicknessOut(double value)     {caseThicknessOut_ = value;}
           
           void FPstyrofoamZLength(double value)          {FPStyrofoamZLength_ = value;}
           void FPCarbonZLength(double value)             {FPCarbonZLength_ = value;}
           void coolFPPipeRadius(double value)            {coolFPPipeRadius_ = value;}
           void nPipes(int value)                         {nPipes_ = value;}
           void pipeRadius(double value)                  {pipeRadius_ = value;}
           void pipeThickness(double value)               {pipeThickness_ = value;}           
	   void pipeTorRadius(std::vector<double>& value) {pipeTorRadius_ = value;}
           void pipeSeparation(double value)              {pipeSeparation_=value;}
          
           void nCrate(int value)                   {nCrate_ = value;}
           void nBoard(int value)                   {nBoard_ = value;}
           void nCrateBeforeSpace(int value)        {nCrateBeforeSpace_ = value;}
           void crateXLength(double value)          {crateXLength_ = value;}
           void crateYLength(double value)          {crateYLength_ = value;}
           void crateZLength(double value)          {crateZLength_ = value;}
           void crateFShieldThickness(double value) {crateFShieldThick_ = value;}
           void crateBShieldThickness(double value) {crateBShieldThick_ = value;}
           void crateTThickness(double value)       {crateTThick_ = value;}
           void crateSThickness(double value)       {crateSThick_ = value;}
           void crateFShieldYLength(double value)   {crateFShieldYLength_ = value;}
           void crateFShieldDeltaZ(double value)    {crateFShieldDeltaZ_ = value;}
           void crateRadiusIn(double value)         {crateRadiusIn_ = value;}
           void cratephi0(double value)             {cratephi0_ = value;}
           void crateDeltaPhi(double value)         {crateDeltaPhi_ = value;}
           void radiatorThickness(double value)     {radiatorThickness_ = value;}
           void activeStripThickness(double value)  {activeStripThickness_ = value;}
           void passiveStripThickness(double value) {passiveStripThickness_ = value;}
           
           
	   
           
           
           
           double envelopeRadiusIn()      const {return envelopeRadiusIn_;}
           double envelopeRadiusOut()     const {return envelopeRadiusOut_;}
           double envelopeZ0()            const {return envelopeZ0_;}
           double envelopeZ1()            const {return envelopeZ1_;}
                      
           bool   crystalShift()          const {return crystalShift_;}
           double crystalZLength()        const {return crystalZLength_;}
           double crystalXYLength()       const {return crystalXYLength_;}
           double crystalVolume()         const {return crystalVolume_;}
           double wrapperThickness()      const {return wrapperThickness_;}
	   double refractiveIndex()       const {return refractiveIndex_; }
	   double crystalDecayTime()      const {return crystalDecayTime_; }
           
	   int    nROPerCrystal()         const {return nROPerCrystal_;}
           double roXLength()             const {return roXLength_;}
           double roYLength()             const {return roYLength_;}
           double roZLength()             const {return roZLength_;}
           double FEEXLength()            const {return FEEXLength_;}
           double FEEYLength()            const {return FEEYLength_;}
           double FEEZLength()            const {return FEEZLength_;}
           double FEEBoxThickness()       const {return FEEBoxThickness_;}
           double BPHoleXLength()         const {return BPHoleXLength_;}
           double BPHoleYLength()         const {return BPHoleYLength_;}
           double BPHoleZLength()         const {return BPHoleZLenght_;}
           double stripThickness()        const {return stripThickness_;}
           double stripYLength()          const {return stripLengthY_;}
           double coolBPPipeRadius()      const {return coolBPPipeRadius_;}

           double outerRingEdgeZLength()  const {return outerRingEdgeZLength_;}
           double outerRingEdgeRLength()  const {return outerRingEdgeRLength_;}
           double caseThicknessIn()       const {return caseThicknessIn_;}
           double caseThicknessOut()      const {return caseThicknessOut_;}
          
           double FPStyrofoamZLength()    const {return FPStyrofoamZLength_;}
           double FPCarbonZLength()       const {return FPCarbonZLength_;}
           double coolFPPipeRadius()      const {return coolFPPipeRadius_;}
           int    nPipes()                const {return nPipes_;}
           double pipeRadius()            const {return pipeRadius_;}
           double pipeThickness()         const {return pipeThickness_;}
           const std::vector<double>& pipeTorRadius() const {return pipeTorRadius_;}
           double pipeSeparation()       const {return pipeSeparation_;}
                
           int    nCrate()                const {return nCrate_;}
           int    nBoard()                const {return nBoard_;}
           int    nCrateBeforeSpace()     const {return nCrateBeforeSpace_;}
           double crateXLength()          const {return crateXLength_;}
	   double crateYLength()          const {return crateYLength_;}
	   double crateZLength()          const {return crateZLength_;}
           double crateFShieldThickness() const {return crateFShieldThick_;}
           double crateBShieldThickness() const {return crateBShieldThick_;}
           double crateTThickness()       const {return crateTThick_;}
           double crateSThickness()       const {return crateSThick_;}
           double crateFShieldYLength()   const {return crateFShieldYLength_;}
           double crateFShieldDeltaZ()    const {return crateFShieldDeltaZ_;}
           double crateRadiusIn()         const {return crateRadiusIn_;}
           double cratephi0()             const {return cratephi0_;}
           double crateDeltaPhi()         const {return crateDeltaPhi_;}          
           double radiatorThickness()     const {return radiatorThickness_;}
           double activeStripThickness()  const {return activeStripThickness_;}
           double passiveStripThickness() const {return passiveStripThickness_;}
           
           
            
           
       private:

          double envelopeRadiusIn_;
          double envelopeRadiusOut_;
          double envelopeZ0_;
          double envelopeZ1_;
          
          bool   crystalShift_;
	  double crystalXYLength_;
	  double crystalZLength_;
	  double crystalVolume_;
          double wrapperThickness_;
          double refractiveIndex_;
          double crystalDecayTime_;

          int    nROPerCrystal_;
          double roXLength_;
          double roYLength_;
          double roZLength_;
          double FEEXLength_;
          double FEEYLength_;
          double FEEZLength_;
          double FEEBoxThickness_;
          double BPHoleXLength_;
          double BPHoleYLength_;
          double BPHoleZLenght_;
          double stripThickness_;
          double stripLengthY_;
          double coolBPPipeRadius_;
          
	  
          double outerRingEdgeZLength_;
          double outerRingEdgeRLength_;
          double caseThicknessIn_;
          double caseThicknessOut_;

          double FPStyrofoamZLength_;
          double FPCarbonZLength_;
          double coolFPPipeRadius_;
	  unsigned int         nPipes_;
	  double               pipeRadius_;
	  double               pipeThickness_;
	  std::vector<double>  pipeTorRadius_;
          double               pipeSeparation_;
          
          int    nCrate_;
          int    nBoard_;
          int    nCrateBeforeSpace_;
          double crateXLength_;
          double crateYLength_;
          double crateZLength_;
          double crateFShieldThick_;
          double crateBShieldThick_;
          double crateTThick_;
          double crateSThick_;
          double crateFShieldYLength_;
          double crateFShieldDeltaZ_;
	  double crateRadiusIn_;
          double cratephi0_;
          double crateDeltaPhi_;
          double radiatorThickness_;
          double activeStripThickness_;
          double passiveStripThickness_;           
     };

}    

#endif 
