#include "TObject.h"
#include "MyVertex.h"



ClassImp(MyVertex)
///////////////////////////////////////////
//  Class containing any kind of vertex  //
///////////////////////////////////////////


//________________________________________________________________________
MyVertex::MyVertex():MyPoint(),
  dmMult(0){
   //default constructor
 }


//___________________________________________________________________________
MyVertex::MyVertex(double X, double Y, double Z, int mult):MyPoint(),
 dmMult(mult){
	//standard constructor 
}	     

//___________________________________________________________________________
MyVertex::MyVertex(const MyVertex& source):MyPoint(source)
{
  dmMult = source.dmMult;
  //copy constructor  
}


//___________________________________________________________________________
MyVertex::~MyVertex()	 {
  //destructor
}

//___________________________________________________________________________
MyVertex& MyVertex::operator=(const MyVertex& source){
  if(this == &source) return *this;
  this->~MyVertex();
  new(this) MyVertex(source);
  return *this;
  //= operator
}