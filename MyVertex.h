#ifndef MyVertex_H
#define MyVertex_H

#include "TObject.h"
#include "MyPoint.h"

class MyVertex : public MyPoint
{

  public:
  
    MyVertex();
    MyVertex(MyPoint* Point, int mult);
    MyVertex(double X, double Y, double Z, int mult);
    MyVertex(const MyVertex& source);
    virtual ~MyVertex();
    MyVertex& operator=(const MyVertex& source);		

    int GetMult() const {return dmMult;}

    void SetMult(int mult) {dmMult = mult;}
  
  private:
    int dmMult;
  

  ClassDef(MyVertex,1)

};


#endif 
