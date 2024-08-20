#include "TObject.h"
#include "TMath.h"
#include "MyPoint.h"



ClassImp(MyPoint)

//________________________________________________________________________
MyPoint::MyPoint():TObject(),
 dmX(0.),
 dmY(0.),
 dmZ(0.){
   // default constructor
 }


//___________________________________________________________________________
MyPoint::MyPoint(double X, double Y, double Z):TObject(),
 dmX(X),
 dmY(Y),
 dmZ(Z){
	//standard constructor 
}	     

//___________________________________________________________________________
MyPoint::MyPoint(const MyPoint& source):TObject(source)
{
  dmX = source.dmX;
  dmY = source.dmY;
  dmZ = source.dmZ;
  //copy constructor  
}

//___________________________________________________________________________
MyPoint::~MyPoint()	 {
  // destructor
}

//___________________________________________________________________________
MyPoint& MyPoint::operator=(const MyPoint& source){
  if(this == &source) return *this;
  this->~MyPoint();
  new(this) MyPoint(source);
  return *this;
  //copy operator
}

//___________________________________________________________________________
double MyPoint::GetRadiusXY() const {
  return TMath::Sqrt(dmX*dmX + dmY*dmY);
}

