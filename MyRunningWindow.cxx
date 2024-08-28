#include "MyRunningWindow.h"

#include "TMath.h"  
#include "Riostream.h"

//using std::vector;

ClassImp(MyRunningWindow)

//___________________________________________________________________________
MyRunningWindow::MyRunningWindow() : TObject(),
 dmSize(0.),
 dmStep(0.){
  //Default constructor
}

//___________________________________________________________________________
MyRunningWindow::MyRunningWindow(double size, double step) : TObject(),
 dmSize(size),
 dmStep(step){
  //Standard constructor
}

//___________________________________________________________________________
MyRunningWindow::MyRunningWindow(const MyRunningWindow& source) : TObject(source)
{
  dmSize = source.dmSize;
  dmStep = source.dmStep;
  //Copy constructor  
}

//___________________________________________________________________________
MyRunningWindow::~MyRunningWindow()
{
  //Default destructor
}

//___________________________________________________________________________
MyRunningWindow& MyRunningWindow::operator=(const MyRunningWindow& source){
    if(this == &source) return *this;
    this->~MyRunningWindow();
    new(this) MyRunningWindow(source);
    return *this;
    //Copy operator
}

//___________________________________________________________________________
double MyRunningWindow::Average(vector<double> vector, double lower_bound, double upper_bound){
  double sum = 0.;
  double counter = 0.;
  int vector_dim = vector.size();
  for(int i = 0; i < vector_dim; i++){
    if(vector.at(i) >= lower_bound && vector.at(i) <= upper_bound){
      sum += vector.at(i);
      counter++;
    }
  }
  return sum/(double)counter;
}

//___________________________________________________________________________
double MyRunningWindow::running_window(vector<double> vector, bool &reconstructed_flag){

  reconstructed_flag = 1; //set to 1 to try to reconstruct the vertex

  int max_counter = 0; //maximum number of "vertexlets" counted in a window
  double reconstructed_vertex = 0.; //reconstructed vertex
  int k_start = 0; //index of the first element of the window
  int i_max = 0; //index of the window with the maximum number of "vertexlets"

  if(vector.empty()){
    reconstructed_flag = 0;
    return 0.;
  }
  else{
    double z0 = vector.at(0) - 0.1; //in cm, starting point of the first window
    int vector_dim = vector.size();
    for(int i = 0; vector.at(vector_dim - 1) > z0 + i*dmStep; i++){
      bool inside_flag = 1; //flag to check if the window is inside the vector
      bool start_flag = 1; //flag to save the index of the first element of the window in the next iteration
      int counter = 0; //counts the number of "vertexlets" in the window
      int k = k_start; //saves the index of the 1st element of the window
      while((k < vector_dim) && (inside_flag)){
        if((vector.at(k) >= z0 + i*dmStep) && (vector.at(k) <= z0 + i*dmStep + dmSize)){
          if(start_flag){
            k_start = k;
            start_flag = 0;
          }
          counter++;
        }
        else if(vector.at(k) > z0 + i*dmStep + dmSize) inside_flag = 0;
        k++;
      }
      if(counter > max_counter){
        max_counter = counter;
        reconstructed_flag = 1;
        i_max = i;
      }
      else if(counter == max_counter && i - i_max != 1){
        reconstructed_flag = 0; //if the maximum is not unique, the vertex is not reconstructed (the second condition means that the windows are not consecutive)
      }
    }

    int left_counts = 0; //counts in the windows before and after the one with the maximum number of "vertexlets"
    int right_counts = 0;
    for(int i = 0; (i < vector_dim) && (vector.at(i) < z0 + (i_max + 1)*dmStep + dmSize); i++) {
      //second condition stops the loop on the vector after the region of interest
      if((vector.at(i) >= z0 + (i_max - 1)*dmStep) && (vector.at(i) <= z0 + i_max*dmStep)) left_counts++;
      if((vector.at(i) >= z0 + i_max*dmStep + dmSize) && (vector.at(i) <= z0 + (i_max + 1)*dmStep + dmSize)) right_counts++;
    }
    
    if (left_counts >= right_counts) reconstructed_vertex = Average(vector, z0 + (i_max - 1)*dmStep, z0 + (i_max)*dmStep + dmSize);
    else if (left_counts < right_counts) reconstructed_vertex = Average(vector, z0 + i_max*dmStep, z0 + (i_max + 1)*dmStep + dmSize);

  }
 
  return reconstructed_vertex;

}