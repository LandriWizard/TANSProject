#ifndef MYPHYSICS_H
#define MYPHYSICS_H

#include <TObject.h>
#include "MyParticle.h"

#include "MyPoint.h"

class MyPhysics : public TObject{

  public:
    MyPhysics();
    MyPhysics(double R, double H);
    MyPhysics(const MyPhysics& source);
    virtual ~MyPhysics();

    double GetR() const {return dmR;}
    double GetH() const {return dmH;}

    void SetR(double r) {dmR = r;}
    void SetH(double h) {dmH = h;}

    MyPoint Hit(MyPoint* Point, MyParticle* Particle);


  private:
    double dmR; //SQRT(X^2+Y^2)
    double dmH; //Lenght of the detectors

  ClassDef(MyPhysics,1)

};
#endif