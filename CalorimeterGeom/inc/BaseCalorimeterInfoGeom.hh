#ifndef CalorimeterGeom_BaseCalorimeterInfoGeom_hh
#define CalorimeterGeom_BaseCalorimeterInfoGeom_hh
//
//
// Contains data for the BaseCalorimeter class
//
// Might consider to have a class for the calibration system and pipes later (or might not...)
//
// Original author B. Echenard
//

#include <vector>

namespace mu2e {

    

    class BaseCalorimeterInfoGeom {


       public:

           BaseCalorimeterInfoGeom(): 
	      _crystalNumEdges(0),_crystalShift(0),_crystalHalfTrans(0),_crystalHalfLength(0),_crystalVolume(0),_wrapperThickness(0),
	      _nROPerCrystal(0),_roHalfTrans(0),_roHalfThickness(0),_roElecHalfX(0),_roElecHalfY(0),_roElecHalfZ(0),
	      _crateRadiusIn(0),_crateRadiusOut(0),_crateHalfLength(0),
	      _caseThickness(0),_envelopeInRadius(0),_envelopeOutRadius(0),_envelopeZ0(0),_envelopeZ1(0),
	      _nPipes(0),_pipeRadius(0),_pipeThickness(0),_pipeTorRadius()	      
	   {}
	     
           virtual ~BaseCalorimeterInfoGeom() {}

           void crystalNedges(int value)           {_crystalNumEdges = value;}
           void crystalShift(bool value)           {_crystalShift = value;}
           void crystalHalfLength(double value)    {_crystalHalfLength = value;}
           void crystalHalfTrans(double value)     {_crystalHalfTrans = value;}
           void crystalVolume(double value)        {_crystalVolume = value;}
	   void wrapperThickness(double value)     {_wrapperThickness = value;}
           
	   void nROPerCrystal(int value)           {_nROPerCrystal = value;}
           void roHalfTrans(double value)          {_roHalfTrans = value;}
           void roHalfThickness(double value)      {_roHalfThickness = value;}
           void roElecHalfX(double value)          {_roElecHalfX = value;}
           void roElecHalfY(double value)          {_roElecHalfY = value;}
           void roElecHalfZ(double value)          {_roElecHalfZ = value;}
           void crateRadiusIn(double value)        {_crateRadiusIn = value;}
           void crateRadiusOut(double value)       {_crateRadiusOut = value;}
           void crateHalfLength(double value)      {_crateHalfLength = value;}
          
	   void caseThickness(double value)        {_caseThickness = value;}
           void envelopeInRadius(double value)     {_envelopeInRadius = value;}
           void envelopeOutRadius(double value)    {_envelopeOutRadius = value;}
           void envelopeZ0(double value)           {_envelopeZ0 = value;}
           void envelopeZ1(double value)           {_envelopeZ1 = value;} 
	   void refractiveIndex(double value)      {_refractiveIndex = value;}
	   void crystalDecayTime(double value)     {_crystalDecayTime = value;}

           
	   
	   int    crystalNedges()       const      {return _crystalNumEdges;}
           bool   crystalShift()        const      {return _crystalShift;}
           double crystalHalfLength()   const      {return _crystalHalfLength;}
           double crystalHalfTrans()    const      {return _crystalHalfTrans;}
           double crystalVolume()       const      {return _crystalVolume;}
           double wrapperThickness()    const      {return _wrapperThickness;}
           
	   int    nROPerCrystal()       const      {return _nROPerCrystal;}
           double roHalfTrans()         const      {return _roHalfTrans;}
           double roHalfThickness()     const      {return _roHalfThickness;}
           double roElecHalfX()         const      {return _roElecHalfX;}
           double roElecHalfY()         const      {return _roElecHalfY;}
           double roElecHalfZ()         const      {return _roElecHalfZ;}
           double crateRadiusIn()       const      {return _crateRadiusIn;}
           double crateRadiusOut()      const      {return _crateRadiusOut;}
           double crateHalfLength()     const      {return _crateHalfLength;}
           
	   double caseThickness()       const      {return _caseThickness;}
           double envelopeInRadius()    const      {return _envelopeInRadius;}
           double envelopeOutRadius()   const      {return _envelopeOutRadius;}
           double envelopeZ0()          const      {return _envelopeZ0;}
           double envelopeZ1()          const      {return _envelopeZ1;}
	   double refractiveIndex()     const      {return _refractiveIndex; }
	   double crystalDecayTime()    const      {return _crystalDecayTime; }


           void nPipes(int value)                  {_nPipes = value;}
           void pipeRadius(double value)           {_pipeRadius = value;}
           void pipeThickness(double value)        {_pipeThickness = value;}           
	   void pipeTorRadius(std::vector<double>& value) { _pipeTorRadius = value;}

           int    nPipes()             const       {return _nPipes;}
           double pipeRadius()         const       {return _pipeRadius;}
           double pipeThickness()      const       {return _pipeThickness;}
           std::vector<double> const& pipeTorRadius() const  {return _pipeTorRadius;}



       private:

          int    _crystalNumEdges;
          bool   _crystalShift;
	  double _crystalHalfTrans;
	  double _crystalHalfLength;
	  double _crystalVolume;
          double _wrapperThickness;

          int    _nROPerCrystal;
          double _roHalfTrans;
          double _roHalfThickness;
          double _roElecHalfX;
          double _roElecHalfY;
          double _roElecHalfZ;
	  
	  double _crateRadiusIn;
	  double _crateRadiusOut;
          double _crateHalfLength;

          double _caseThickness;
          double _envelopeInRadius;
          double _envelopeOutRadius;
          double _envelopeZ0;
          double _envelopeZ1;
	   
          double _refractiveIndex;
          double _crystalDecayTime;

	  unsigned int         _nPipes;
	  double               _pipeRadius;
	  double               _pipeThickness;
	  std::vector<double>  _pipeTorRadius;


     };

}    

#endif 
