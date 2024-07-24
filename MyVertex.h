#ifndef MyVertex_H
#define MyVertex_H

#include "MyPoint.h"

class MyVertex : public MyPoint
{

  public:
  
    MyVertex();
    MyVertex(MyPoint* Point, int mult);
//Inserisco MyVertex(double X, double Y, double Z, int mult);   ??
    MyVertex(const MyVertex& source);
    virtual ~MyVertex();
    MyVertex& operator=(const MyVertex& source);		

//GETTERS  
    int GetMult() const {return dmMult;}

//SETTERS
    void SetMult(int mult) {dmMult = mult;}
  
  private:
  
    int dmMult;
  
  ClassDef(MyVertex,1)

};


#endif 
