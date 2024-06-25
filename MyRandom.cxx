#include <Riostream.h>
#include <TMath.h>
#include "MyRandom.h"

ClassImp(MyRandom)
///////////////////////////////////////////
// Class used to generate random numbers //
///////////////////////////////////////////

//-----------------------------------------------------------//
MyRandom::MyRandom(): TRandom3(),
fAlpha(0.),
fBig(0.),
fPi(0.) {
  //Default constructor
}

//-----------------------------------------------------------//
MyRandom::MyRandom(double alpha, unsigned int seed): TRandom3(seed),
  fAlpha(alpha),
  fBig(0.),
  fPi(0.) {
  //Standard constructor
  Initialisation();
}

//-----------------------------------------------------------//
MyRandom::~MyRandom()
{
  //Default destructor
}

//-----------------------------------------------------------//
double MyRandom::Initialisation()
{
  fPi = TMath::Pi();     //Definition of Pi, used in many places in the class
  fBig = 1./fAlpha;      //Maximum for the function implemented in Func()
}

//-----------------------------------------------------------//
double MyRandom::Func(double x)
{
  return 1./(TMath::Sin(x)*TMath::Sin(x) + fAlpha*TMath::Cos(x)*TMath::Cos(x));
}
 
//-----------------------------------------------------------//
double MyRandom::Rejection()
{
  double x, y;
  do{
    x = 2.*Rndm()*TMath::Pi();
    y = Rndm()*fBig;
  }while(y > Func(x));
  return x;
}
