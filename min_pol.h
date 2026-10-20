#ifndef MIN_POL_H
#define MIN_POL_H
#include<iostream>
class Min_Pol{
  public:
  
  //Constructor & Descructor

  Min_Pol();
  Min_Pol(int size);
  ~Min_Pol();
  double& operator[](const int &i);
  const double& operator[](const int &i) const;
  Min_Pol& operator=(const Min_Pol& rhs);
  int get_size() const; 
  private:
  
  double* coeffecients;
  int size;
};
Min_Pol::Min_Pol(){
  size = 0;
  coeffecients = nullptr;
}
Min_Pol::Min_Pol(int size){
  coeffecients = new double[size];
  for(int i = 0 ; i < size; ++i){
    coeffecients[i] = 0.0;
  }
  this->size = size;
}

Min_Pol::~Min_Pol(){

  delete[] coeffecients;
}
double& Min_Pol::operator[](const int &i){
  if(i < 0 || i >= size){
   throw std::runtime_error("i out of range for Min_Pol"); 
  }

  return coeffecients[i];
}
const double& Min_Pol::operator[](const int &i) const{
   if(i < 0 || i >= size){
   throw std::runtime_error("i out of range for Min_Pol"); 
  }

  return coeffecients[i]; 
}
Min_Pol& Min_Pol::operator=(const Min_Pol& rhs){
  if(this == &rhs){
    return *this;
  }

  delete[] coeffecients;

  size = rhs.size;

  coeffecients = new double[size];
  for(int i = 0; i < size; ++i){
    coeffecients[i] = rhs[i];
  }
  return *this;
}

int Min_Pol::get_size() const{
  return size;
}


void print_minimal_polynomial( const Min_Pol &mp){
  int size = mp.get_size();
  printf("A^%d = ", size);
  for(int i = size - 1; i >=0; --i){
    printf("%f A^%d", mp[i], i);
    if( i != 0){
      printf(" + ");
    }
  }
  printf("\n");
}
#endif
