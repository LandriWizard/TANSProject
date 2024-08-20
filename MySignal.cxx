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
MySignal::MySignal(double r, double z, double Phi, int particle_flag): TObject(),
  dmR(r),
  dmZ(z),
  dmPhi(Phi),
  dmParticleFlag(particle_flag){
  //Standard constructor
}

MySignal::MySignal(MyPoint* Point, int particle_flag): TObject()
{

  double x = Point->GetX();
  double y = Point->GetY();
  double r = Point->GetRadiusXY();
  double tmp;

  if(y >= 0.) tmp = TMath::ACos(x/r);
  else tmp = 2.*TMath::Pi() - TMath::ACos(x/r);

  dmR = r;
  dmZ = Point->GetZ();
  dmPhi = tmp;
  dmParticleFlag = particle_flag;
}

//___________________________________________________________________________
MySignal::MySignal(const MySignal& source): TObject(source)
{
  dmR = source.dmR;
  dmZ = source.dmZ;
  dmPhi = source.dmPhi;
  dmParticleFlag = source.dmParticleFlag;
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
