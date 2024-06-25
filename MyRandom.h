#ifndef MYRANDOM_H
#define MYRANDOM_H

#include <TRandom3.h>

class MyRandom : public TRandom3 {
//////////////////////////////////////////
//  Class for random number generators  //
//////////////////////////////////////////

  public:
    MyRandom();
    MyRandom(double alpha, unsigned int seed);
    virtual ~MyRandom();
    double Initialisation();
    double Func(double x);
    double Rejection();





  private:
    double fAlpha; //parameter in the functions
    double fBig;   //upper limit for the rejection method
    double fPi;    //Pi, self explainatory



};
#endif
