#ifndef MYPHYSICS_H
#define MYPHYSICS_H

#include "TObject.h"

#include "MySignal.h"
#include "MyParticle.h"
#include "MyPoint.h"
#include "MyRandom.h"


class MyPhysics : public TObject{

  public:
    MyPhysics();
    MyPhysics(double R, double H, double scattering_theta, double smearing_z, double smearing_rphi);
    MyPhysics(const MyPhysics& source);
    virtual ~MyPhysics();
    MyPhysics& operator=(const MyPhysics& source);

    double GetR() const {return dmR;}
    double GetH() const {return dmH;}
    double GetScatteringAngle() const {return dmScatteringTheta;}

    void SetR(double r) {dmR = r;}
    void SetH(double h) {dmH = h;}
    void SetScatteringAngle(double theta) {dmScatteringTheta = theta;}

    MyPoint Transport(MyPoint* Point, MyParticle* Particle);

    MyParticle NoScattering(MyParticle* Particle) {return *Particle;}
    MyParticle MultipleScattering(MyParticle* Particle);

    MySignal NoSmearing(MySignal* Signal) {return *Signal;}
    MySignal SmearingOn(MySignal* Signal);



  private:
    double dmR; //SQRT(X^2+Y^2)
    double dmH; //Lenght of the detectors
    double dmScatteringTheta; //scattering angle rms
    double dmSmearingZ;
    double dmSmearingRPhi;

  ClassDef(MyPhysics,1)

};
#endif