#include "MyPhysics.h"

#include "TMath.h"
#include "Riostream.h"

ClassImp(MyPhysics)

////////////////////////////////////////////////
// Class containing all the physics phenomena //
// that we need to simulate                   //
////////////////////////////////////////////////

//___________________________________________________________________________
MyPhysics::MyPhysics(): TObject(),
  dmR(0.),
  dmH(0.){
  //Default constructor
}

//___________________________________________________________________________
MyPhysics::MyPhysics(double R, double H): TObject(),
  dmR(R),
  dmH(H){
  //Standard constructor
}

//___________________________________________________________________________
MyPhysics::MyPhysics(const MyPhysics& source): TObject(source)
{
  dmR = source.dmR;
  dmH = source.dmH;
  //copy constructor  
}

//___________________________________________________________________________
MyPhysics::~MyPhysics()
{
  //Default destructor
}

//___________________________________________________________________________
MyPoint MyPhysics::Hit(MyPoint* Point, MyParticle* Particle){
  double x0 = Point->GetX();
  double y0 = Point->GetY();
  double z0 = Point->GetZ();

  double c1 = TMath::Sin(Particle->GetTheta())*TMath::Cos(Particle->GetPhi());
  double c2 = TMath::Sin(Particle->GetTheta())*TMath::Sin(Particle->GetPhi());
  double c3 = TMath::Cos(Particle->GetTheta());

  double delta = (x0*c1 + y0*c2)*(x0*c1 + y0*c2) - (c1*c1 + c2*c2)*(x0*x0* + y0*y0 - dmR*dmR);

  double tminus = (-x0*c1 - y0*c2 - TMath::Sqrt(delta))/(c1*c1 + c2*c2);
  double tplus = (-x0*c1 - y0*c2 + TMath::Sqrt(delta))/(c1*c1 + c2*c2);

  double t;
  if(tplus > 0.) t = tplus;
  else t = tminus;

  std::cout << "t_minus = " << tminus << ", t_plus = " << tplus << ", t = " << t << ";" << std::endl; 

  double x = x0 + c1*t;
  double y = y0 + c2*t;
  double z = z0 + c3*t;

  MyPoint hit(x,y,z);

  return hit;
}
