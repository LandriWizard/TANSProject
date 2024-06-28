#include "TObject.h"
#include "MyVertex.h"



ClassImp(MyVertex)

//________________________________________________________________________
MyVertex::MyVertex():MyPoint(),
 dmPoint(NULL),
 dmMult(0){
   // default constructor
 }


//___________________________________________________________________________
MyVertex::MyVertex(MyPoint* Point, int mult):MyPoint(),
 dmPoint(Point),
 dmMult(mult){
	//standard constructor 
}	     

//___________________________________________________________________________
MyVertex::MyVertex(const MyVertex& source):MyPoint(source)
{
  dmPoint = source.dmPoint;
  dmMult = source.dmMult;
  //copy constructor  
}


//___________________________________________________________________________
MyVertex::~MyVertex()	 {
  // destructor
}

//___________________________________________________________________________
MyVertex& MyVertex::operator=(const MyVertex& source){
  if(this == &source) return *this;
  this->~MyVertex();
  new(this) MyVertex(source);
  return *this;
  //= operator
}