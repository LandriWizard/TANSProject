#ifndef MYRANDOM_H
#define MYRANDOM_H

#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>

#include "MyVertex.h"

class MyRandom : public TRandom3 {
//////////////////////////////////////////
//  Class for random number generators  //
//////////////////////////////////////////

  public:
    MyRandom();
    MyRandom(const char* input_file, unsigned int seed);
    MyRandom(const MyRandom& source);
    virtual ~MyRandom();
    MyRandom& operator=(const MyRandom& source);

    int RndmMult();
    double RndmTheta();

    static bool GetFlag() {return dmFileFlag;}

  private:
    static bool dmFileFlag; //Flag for the existence of input_file; TRUE if input_file does not exist

    TH1D* dmMult;
    TH1D* dmEta;

};
#endif
