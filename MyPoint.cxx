#include "TObject.h"
#include "MyPoint.h"



ClassImp(MyPoint)

//________________________________________________________________________
MyPoint::MyPoint():TObject(),
 fX(0.),
 fY(0.),
 fZ(0.){
   // default constructor
 }


//___________________________________________________________________________
MyPoint::MyPoint(double X, double Y, double Z):TObject(),
 fX(X),
 fY(Y),
 fZ(Z){
	//standard constructor 
}	     

//___________________________________________________________________________
MyPoint::~MyPoint()	 {
  // destructor
}
