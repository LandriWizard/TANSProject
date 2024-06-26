#include "TObject.h"
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
