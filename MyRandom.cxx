#include <Riostream.h>
#include <TMath.h>
#include "MyRandom.h"

ClassImp(MyRandom)
///////////////////////////////////////////
// Class used to generate random numbers //
///////////////////////////////////////////

//-----------------------------------------------------------//
MyRandom::MyRandom(): TRandom3(),
  dmMult(0){
  //Default constructor
}

//-----------------------------------------------------------//
MyRandom::MyRandom(const char* input_file,unsigned int seed): TRandom3(seed),
  dmMult(){
  //Standard constructor
}

//___________________________________________________________________________
MyRandom::MyRandom(const MyRandom& source):TRandom3(source)
{
  dmMult = source.dmMult;
  //copy constructor  
}


//-----------------------------------------------------------//
MyRandom::~MyRandom()
{
  //Default destructor
}

//___________________________________________________________________________
MyRandom& MyRandom::operator=(const MyRandom& source){
  if(this == &source) return *this;
  this->~MyRandom();
  new(this) MyRandom(source);
  return *this;
  //= operator
}