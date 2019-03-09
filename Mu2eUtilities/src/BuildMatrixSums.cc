//S Middleton - March 19
//Class to store the summations which for elements of gamma and beta matrices for straight track cosmic fits
#include <cstdio>
#include <math.h>

#include "DataProducts/inc/XYZVec.hh"
#include "Mu2eUtilities/inc/BuildMatrixSums.hh"

BuildMatrixSums::BuildMatrixSums() {
  clear();
}

BuildMatrixSums::BuildMatrixSums(const BuildMatrixSums& S) {
  init(S);
 
}

BuildMatrixSums::~BuildMatrixSums() {
}
void BuildMatrixSums::clear(){
  
  betaX00 = 0;
  betaX10 = 0;
  gammaX00 = 0;
  gammaX01 = 0;
  gammaX11 = 0;
  betaY00 =  0;
  betaY10 =  0;
  gammaY00 = 0;
  gammaY01 = 0;
  gammaY11 = 0;
 
}
void BuildMatrixSums::init(const BuildMatrixSums& S) {
  
  betaX00 = S.betaX00;
  betaX10 = S.betaX10;
  gammaX00 = S.gammaX00;
  gammaX01 = S.gammaX01;
  gammaX11 = S.gammaX11;
  betaY00 =  S.betaY00;
  betaY10 =  S.betaY10;
  gammaY00 = S.gammaY00;
  gammaY01 = S.gammaY01;
  gammaY11 = S.gammaY11;
 
}

void BuildMatrixSums::addPoint(int i, XYZVec point_i, XYZVec track_direction, double errX, double errY){
       
        XYZVec track_y(0,track_direction.y(),0);
        XYZVec track_x(track_direction.x(),0,0);
	//For BetaX:
	betaX00 += point_i.Dot(track_x)/pow(errX,2);
        betaX10 += point_i.z()*point_i.Dot(track_x)/pow(errX,2);
	//For GammaX:
	gammaX00 +=pow(errX,2);
        gammaX01  += point_i.z()/pow(errX,2);
	gammaX11 += pow(point_i.z(),2) /pow(errX,2);
	//For BetaY:
	betaY00 += point_i.Dot(track_y)/pow(errY,2);
        betaY10 += point_i.z()*point_i.Dot(track_y)/pow(errY,2);
	//For GammaY:
	gammaY00 +=pow(errY,2);
        gammaY01  += point_i.z()/pow(errY,2);
	gammaY11 += pow(point_i.z(),2) /pow(errY,2);
  
}

void BuildMatrixSums::removePoint(XYZVec point_i, XYZVec track_direction, double errX, double errY){
        XYZVec track_y(0,track_direction.y(),0);
        XYZVec track_x(track_direction.x(),0,0);
	//For BetaX:
	betaX00 -= point_i.Dot(track_x)/pow(errX,2);
        betaX10 -= point_i.z()*point_i.Dot(track_x)/pow(errX,2);
	//For GammaX:
	gammaX00 -=pow(errX,2);
        gammaX01  -= point_i.z()/pow(errX,2);
	gammaX11 -= pow(point_i.z(),2) /pow(errX,2);
	//For BetaY:
	betaY00 -= point_i.Dot(track_y)/pow(errY,2);
        betaY10 -= point_i.z()*point_i.Dot(track_y)/pow(errY,2);
	//For GammaY:
	gammaY00 -=pow(errY,2);
        gammaY01  -= point_i.z()/pow(errY,2);
	gammaY11 -= pow(point_i.z(),2) /pow(errY,2);
  
}


double BuildMatrixSums::Get2DParameter(int i, TMatrixD Alpha){
        
	return Alpha[i][0];

}


TMatrixD BuildMatrixSums::GetGammaX(){
	TMatrixD Gamma(2,2);
	Gamma[0][0] = gammaX00;
	Gamma[1][0] = gammaX01;
	Gamma[0][1] = gammaX01;
	Gamma[1][1] = gammaX11;
	return Gamma;
}

TMatrixD BuildMatrixSums::GetBetaX(){
	TMatrixD Beta(2,1);
	Beta[0][0] = betaX00;
	Beta[1][0] = betaX10;
	return Beta;
}

TMatrixD BuildMatrixSums::GetAlphaX(){
        TMatrixD Gamma = GetGammaX();
	TMatrixD Beta  = GetBetaX();
	double* det = NULL;                  
	Gamma.Invert(det);
        TMatrixD Alpha(Gamma*Beta);
	return Alpha;
}

TMatrixD BuildMatrixSums::GetGammaY(){
	TMatrixD Gamma(2,2);
	Gamma[0][0] = gammaY00;
	Gamma[1][0] = gammaY01;
	Gamma[0][1] = gammaY01;
	Gamma[1][1] = gammaY11;
	return Gamma;
}

TMatrixD BuildMatrixSums::GetBetaY(){
	TMatrixD Beta(2,1);
	Beta[0][0] = betaY00;
	Beta[1][0] = betaY10;
	return Beta;
}

TMatrixD BuildMatrixSums::GetAlphaY(){
        TMatrixD Gamma = GetGammaY();
	TMatrixD Beta  = GetBetaY();
	double* det = NULL;                  
	Gamma.Invert(det);
        TMatrixD Alpha(Gamma*Beta);
	return Alpha;
}
