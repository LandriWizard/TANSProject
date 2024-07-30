#ifndef MYSIGNAL_H
#define MYSIGNAL_H

#include "MyPoint.h"

class MySignal : public MyPoint{

  public:
    MySignal();
    MySignal(double r, double z, double Phi);
    MySignal(MyPoint *Point);
    MySignal(const MySignal& source);
    virtual ~MySignal();
    MySignal& operator=(const MySignal& source);

    double GetR() const {return dmR;}
    double GetPhi() const {return dmPhi;}

    void SetR(double r) {dmR = r;}
    void SetH(double Phi) {dmPhi = Phi;}



  private:
    double dmR;
    double dmPhi; //Azimuth angle of the point

  ClassDef(MySignal,1)

};
#endif