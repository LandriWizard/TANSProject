#ifndef MYSIGNAL_H
#define MYSIGNAL_H

#include <TObject.h>

#include "MyPoint.h"

class MySignal : public TObject{

  public:
    MySignal();
    MySignal(double z, double Phi);
    MySignal(MyPoint *Point);
    MySignal(const MySignal& source);
    virtual ~MySignal();

    double GetZ() const {return dmZ;}
    double GetPhi() const {return dmPhi;}

    void SetR(double z) {dmZ = z;}
    void SetH(double Phi) {dmPhi = Phi;}



  private:
    double dmZ; 
    double dmPhi; //Azimuth angle of the point

  ClassDef(MySignal,1)

};
#endif