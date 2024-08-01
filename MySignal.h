#ifndef MYSIGNAL_H
#define MYSIGNAL_H

#include "TObject.h"
#include "MyPoint.h"

class MySignal : public TObject{

  public:
    MySignal();
    MySignal(double r, double z, double Phi, int particle_flag);
    MySignal(MyPoint *Point, int particle_flag);
    MySignal(const MySignal& source);
    virtual ~MySignal();
    MySignal& operator=(const MySignal& source);

    double GetR() const {return dmR;}
    double GetZ() const {return dmZ;}
    double GetPhi() const {return dmPhi;}
    int GetFlag() const {return dmParticleFlag;}

    void SetR(double r) {dmR = r;}
    void SetZ(double z) {dmZ = z;}
    void SetH(double Phi) {dmPhi = Phi;}
    void SetFlag(int particle_flag) {dmParticleFlag = particle_flag;}


  private:
    double dmR;
    double dmZ;
    double dmPhi; //Azimuth angle of the point
    int dmParticleFlag; //Flag used to 

  ClassDef(MySignal,1)

};
#endif