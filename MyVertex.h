#ifndef MyVertex_H
#define MyVertex_H

#include "MyPoint.h"

class MyVertex : public MyPoint
{

  public:
  
    MyVertex();
    MyVertex(MyPoint* Point, int mult);
    MyVertex(const MyVertex& source);
    virtual ~MyVertex();
    MyVertex& operator=(const MyVertex& source);		

//GETTERS  
    MyPoint *GetPoint() const {return dmPoint;}
    int GetMult() const {return dmMult;}

//SETTERS
    void SetPoint(MyPoint *Point) {dmPoint = Point;}
    void SetMult(int mult) {dmMult = mult;}
  
  private:
  
    MyPoint* dmPoint;
    int dmMult;
  
  ClassDef(MyVertex,1)
};


#endif 
