#include "MyRunningWindow.h"

#include "TMath.h"  
#include "Riostream.h"

ClassImp(MyRunningWindow)

///////////////////////////////////
// Class with the implementation //
// of the running window method  //
// for the vertex reconstruction //
///////////////////////////////////

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
  for(int i = 0; i < vector_dim; i++){  //scans the whole vector
    if(vector.at(i) >= lower_bound && vector.at(i) <= upper_bound){ //takes only elements in the wanted interval
      sum += vector.at(i);  //sums the elements in the interval
      counter++;  //counts the elements in the interval
    }
  }
  return sum/(double)counter; //returns the average of the elements in the interval
}

//___________________________________________________________________________
double MyRunningWindow::running_window(vector<double> vector, bool &reconstructed_flag){

  reconstructed_flag = 1; //set to 1 to try to reconstruct the vertex

  int max_counter = 0; //maximum number of "vertexlets" counted in a window
  double reconstructed_vertex = 0.; //reconstructed vertex
  int k_start = 0; //index of the first element of the window
  int i_max = 0; //index of the window with the maximum number of "vertexlets"

  if(vector.empty()){   //stops the function if the vector is empty
    reconstructed_flag = 0;
    return 0.;
  }
  else{
    double z0 = vector.at(0) - 0.1; //in cm, starting point of the first window
    int vector_dim = vector.size();
    for(int i = 0; vector.at(vector_dim - 1) > z0 + i*dmStep; i++){ //scanning the whole vector window by window
      bool inside_flag = 1; //flag to check if the window is inside the vector
      bool start_flag = 1; //flag to save the index of the first element of the window in the next iteration
      int counter = 0; //counts the number of "vertexlets" in the window
      int k = k_start; //saves the index of the 1st element of the window
      while((k < vector_dim) && (inside_flag)){
        if((vector.at(k) >= z0 + i*dmStep) && (vector.at(k) <= z0 + i*dmStep + dmSize)){  //checks if we still are in the window
          if(start_flag){
            k_start = k;    //saves the index of the first element of the window
            start_flag = 0; //makes sure that only the first element is saved
          }
          counter++;  //counting the "vertexlets" in the window
        }
        else if(vector.at(k) > z0 + i*dmStep + dmSize) inside_flag = 0; //if we are out of the window, the loop stops
        k++;  //goes to the next element of the vector
      }
      if(counter > max_counter){  //checks if the number of "vertexlets" in the window is the maximum
        max_counter = counter;    //saves the maximum number of "vertexlets"
        reconstructed_flag = 1;   //if the maximum is reached, the vertex is reconstructed
        i_max = i;                //saves the index of the window with the maximum number of "vertexlets"
      }
      else if(counter == max_counter && i - i_max != 1){
        reconstructed_flag = 0; //if the maximum is not unique, the vertex is not reconstructed (the second condition means that the windows are not consecutive)
      }
    }

    int left_counts = 0; //counts in the windows before and after the one with the maximum number of "vertexlets"
    int right_counts = 0; 
    for(int j = 0; (j < vector_dim) && (vector.at(j) < z0 + (i_max + 1)*dmStep + dmSize); j++) {
      //first condition checks if we are in the vector
      //second condition stops the loop on the vector after the region of interest

      //counts the "vertexlets" in the window before the one with the maximum number of "vertexlets"
      if((vector.at(j) >= z0 + (i_max - 1)*dmStep) && (vector.at(j) <= z0 + i_max*dmStep)) left_counts++; 
      //counts the "vertexlets" in the window after the one with the maximum number of "vertexlets"
      if((vector.at(j) >= z0 + i_max*dmStep + dmSize) && (vector.at(j) <= z0 + (i_max + 1)*dmStep + dmSize)) right_counts++;
    }
    
    //reconstructs the vertex based on the number of "vertexlets" in the windows before and after the one with the maximum number of "vertexlets"
    //this is done by averaging the "vertexlets" in the window with the maximum number of "vertexlets" and the window with the maximum number of "vertexlets" +/- 1
    if (left_counts >= right_counts) reconstructed_vertex = Average(vector, z0 + (i_max - 1)*dmStep, z0 + (i_max)*dmStep + dmSize);
    else if (left_counts < right_counts) reconstructed_vertex = Average(vector, z0 + i_max*dmStep, z0 + (i_max + 1)*dmStep + dmSize);

  }
 
  return reconstructed_vertex;

}