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

    MyPoint Hit(MyPoint* Point, MyParticle* Particle);

  private:
    double dmR; //SQRT(X^2+Y^2)
    double dmH; //Lenght of the detectors

  ClassDef(MyPhysics,1)

};
#endif