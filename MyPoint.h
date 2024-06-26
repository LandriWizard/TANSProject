#ifndef MyPoint_H
#define MyPoint_H

#include "TObject.h"

class MyPoint : public TObject
{

public:

MyPoint();
MyPoint(double X, double Y, double Z);

virtual ~MyPoint();

double GetX() const {return dmX;}; 
double GetY() const {return dmY;};
double GetZ() const {return dmZ;};

double SetX(double x) {dmX = x;};
double SetY(double y) {dmY = y;};
double SetZ(double z) {dmZ = z;};


private:


double dmX;
double dmY;
double dmZ;

ClassDef(MyPoint,1)
};


#endif 
