#ifndef MyPoint_H
#define MyPoint_H

#include "TObject.h"

class MyPoint : public TObject
{

public:

MyPoint();
MyPoint(double X, double Y, double Z);

virtual ~MyPoint();

 double GetX() const {return fX;} 
 double GetY() const {return fY;}
 double GetZ() const {return fZ;}


private:


double fX;
double fY;
double fZ;

ClassDef(MyPoint,1)
};


#endif 
