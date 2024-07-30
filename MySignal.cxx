#include "MySignal.h"

#include "TMath.h"
#include "Riostream.h"

ClassImp(MySignal)

////////////////////////////////////////////////
// Class containing the signals, i.e. the     //
// recorded hits used for the reconstruction  //
////////////////////////////////////////////////

//___________________________________________________________________________
MySignal::MySignal(): TObject(),
  dmR(0.),
  dmZ(0.),
  dmPhi(0.){
  //Default constructor
}

//___________________________________________________________________________
MySignal::MySignal(double r, double z, double Phi): TObject(),
  dmR(r),
  dmZ(z),
  dmPhi(Phi){
  //Standard constructor
}

MySignal::MySignal(MyPoint* Point): TObject()
{

  double x = Point->GetX();
  double y = Point->GetY();
  double r = Point->GetRadius();
  double tmp;

  if(y >= 0.) tmp = TMath::ACos(x/r);
  else tmp = 2.*TMath::Pi() - TMath::ACos(x/r);

  dmR = r;
  dmZ = Point->GetZ();
  dmPhi = tmp;
}

//___________________________________________________________________________
MySignal::MySignal(const MySignal& source): TObject(source)
{
  dmR = source.dmR;
  dmZ = source.dmZ;
  dmPhi = source.dmPhi;
  //copy constructor  
}

//___________________________________________________________________________
MySignal::~MySignal()
{
  //Default destructor
}

//___________________________________________________________________________
MySignal& MySignal::operator=(const MySignal& source){
    if(this == &source) return *this;
    this->~MySignal();
    new(this) MySignal(source);
    return *this;
    //copy operator
}
