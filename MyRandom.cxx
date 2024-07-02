#include <Riostream.h>
#include <TMath.h>
#include "MyRandom.h"

ClassImp(MyRandom)
///////////////////////////////////////////
// Class used to generate random numbers //
///////////////////////////////////////////

bool MyRandom::dmFileFlag = false;

//___________________________________________________________________________
MyRandom::MyRandom(): TRandom3(),
  dmMult(NULL),
  dmEta(NULL){
  //Default constructor
}

//___________________________________________________________________________
MyRandom::MyRandom(const char* input_file,unsigned int seed): TRandom3(seed)
{
  TFile infile(input_file);
  if(infile.IsZombie()) dmFileFlag = true;
  else{
    TH1D* tempEta = (TH1D*)infile.Get("heta2");
    TH1D* tempMult = (TH1D*)infile.Get("hm");
    tempEta->SetDirectory(0);
    tempMult->SetDirectory(0);
    infile.Close();
    dmEta = tempEta;
    dmMult = tempMult;
  }
  infile.Close();
  //Standard constructor
}

//___________________________________________________________________________
MyRandom::MyRandom(const MyRandom& source):TRandom3(source)
{
  dmMult = source.dmMult;
  dmEta = source.dmEta;
  //copy constructor  
}


//___________________________________________________________________________
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

int MyRandom::RndmMult() {
  return dmMult->GetRandom();
}

double MyRandom::RndmTheta(){
  double eta = dmEta->GetRandom();
  return 2.*TMath::ATan(TMath::Exp(-eta));
}
