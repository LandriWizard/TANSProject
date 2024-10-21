#include <Riostream.h>
#include <TMath.h>

#include "MyParticle.h"

ClassImp(MyParticle)

///////////////////////////////////////////////////
// Class containing the "particles" i.e. objects //
// defined by a couple of angles theta, phi that //
// are the direction of the trajectory of the    //
// particle.                                     //
///////////////////////////////////////////////////


//________________________________________________________________________
MyParticle::MyParticle():TObject(),
  dmTheta(0.),
  dmPhi(0.){
    // default constructor
}


//___________________________________________________________________________
MyParticle::MyParticle(double Theta, double Phi):TObject(),
  dmTheta(Theta),
  dmPhi(Phi){
 	//standard constructor 
}	     

//___________________________________________________________________________
MyParticle::MyParticle(const MyParticle& source):TObject(source)
{
  dmTheta = source.dmTheta;
  dmPhi = source.dmPhi;
  //copy constructor  
}

//___________________________________________________________________________
MyParticle::~MyParticle()	 {
  // destructor
}

//___________________________________________________________________________
MyParticle& MyParticle::operator=(const MyParticle& source){
  if(this == &source) return *this;
  this->~MyParticle();
  new(this) MyParticle(source);
  return *this;
  //Copy operator
}