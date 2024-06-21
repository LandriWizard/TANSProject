#ifndef PUNTO_H
#define PUNTO_H

#include "TObject.h"

class Punto : public TObject
{

public:

Punto();
Punto(double X, double Y, double Z);

virtual ~Punto();

 double GetX() const {return fX;} 
 double GetY() const {return fY;}
 double GetZ() const {return fZ;}


private:


double fX;
double fY;
double fZ;

ClassDef(Punto,1)
};


#endif 
