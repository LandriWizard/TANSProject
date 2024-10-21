#include "MyTracklet.h"

ClassImp(MyTracklet)

//////////////////////////////////////////////
// Class containing the Tracklets, i.e. the //
// tracks obtained from hits on layer 2 and //
// compatible hits on layer 1               //
//////////////////////////////////////////////

//___________________________________________________________________________
MyTracklet::MyTracklet(): TObject(),
  dmR1(0.),
  dmR2(0.),
  dmZ1(0.),
  dmZ2{0.}{
  //Default constructor
}

//___________________________________________________________________________
MyTracklet::MyTracklet(double R1, double R2, double Z1, double Z2): TObject(),
dmR1(R1),
dmR2(R2),
dmZ1(Z1),
dmZ2{Z2}{
  //Standard constructor
}

//___________________________________________________________________________
MyTracklet::MyTracklet(MySignal* InnerSignal, MySignal* OuterSignal){
  dmR1 = InnerSignal->GetR();
  dmR2 = OuterSignal->GetR();
  dmZ1 = InnerSignal->GetZ();
  dmZ2 = OuterSignal->GetZ();
  //Alternative constructor using MySignal objects
}


//___________________________________________________________________________
MyTracklet::MyTracklet(const MyTracklet& source): TObject(source)
{
  dmR1 = source.dmR1;
  dmR2 = source.dmR2;
  dmZ1 = source.dmZ1;
  dmZ2 = source.dmZ2;
  //copy constructor  
}

//___________________________________________________________________________
MyTracklet::~MyTracklet()
{
  //Default destructor
}

//___________________________________________________________________________
MyTracklet& MyTracklet::operator=(const MyTracklet& source){
    if(this == &source) return *this;
    this->~MyTracklet();
    new(this) MyTracklet(source);
    return *this;
    //copy operator
}

//___________________________________________________________________________
//Function used to find the intersection between tracklet and the z axis
double MyTracklet::Intersection(){
  double m = (dmR2 - dmR1)/(dmZ2 - dmZ1);
  double z = -(dmR1 - m*dmZ1)/m;
  return z;
}