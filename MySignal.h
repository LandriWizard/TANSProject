#ifndef MYSIGNAL_H
#define MYSIGNAL_H

#include "TObject.h"
#include "MyPoint.h"

class MySignal : public TObject{

  public:
    MySignal();
    MySignal(double r, double z, double Phi);
    MySignal(MyPoint *Point);
    MySignal(const MySignal& source);
    virtual ~MySignal();
    MySignal& operator=(const MySignal& source);

    double GetR() const {return dmR;}
    double GetZ() const {return dmZ;}
    double GetPhi() const {return dmPhi;}

    void SetR(double r) {dmR = r;}
    void SetZ(double z) {dmZ = z;}
    void SetH(double Phi) {dmPhi = Phi;}



  private:
    double dmR;
    double dmZ;
    double dmPhi; //Azimuth angle of the point

  ClassDef(MySignal,1)

};
#endif