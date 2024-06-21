#include "TObject.h"
#include "Punto.h"



ClassImp(Punto)

//________________________________________________________________________
Punto::Punto():TObject(),
 fX(0.),
 fY(0.),
 fZ(0.){
   // default constructor
 }


//___________________________________________________________________________
Punto::Punto(double X, double Y, double Z):TObject(),
 fX(X),
 fY(Y),
 fZ(Z){
	//standard constructor 
}	     

//___________________________________________________________________________
Punto::~Punto()	 {
  // destructor
}
