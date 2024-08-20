#ifndef MyPoint_H
#define MyPoint_H

#include "TObject.h"

class MyPoint : public TObject
{

  public:

    MyPoint();
    MyPoint(double X, double Y, double Z);
    MyPoint(const MyPoint& source);
    virtual ~MyPoint();
    MyPoint& operator=(const MyPoint& source);		

    double GetX() const {return dmX;} 
    double GetY() const {return dmY;}
    double GetZ() const {return dmZ;}
    double GetRadiusXY() const;

    void SetX(double x) {dmX = x;}
    void SetY(double y) {dmY = y;}
    void SetZ(double z) {dmZ = z;}


  private:

    double dmX;
    double dmY;
    double dmZ;

  ClassDef(MyPoint,1)
};


#endif 
