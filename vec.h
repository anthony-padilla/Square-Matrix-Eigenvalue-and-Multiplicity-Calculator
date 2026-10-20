#ifndef VEC_H
#define VEC_H
#include <iostream>
#include <cmath>
class vec{

  public:
    ~vec();
    vec(int i);

    int get_size() const;
    double& operator[](int i);
    const double& operator[](int i)const;
    vec operator+(const vec& rhs);
    vec operator-(const vec& rhs);
    vec operator*(const double s);
    vec& operator=(const vec& other);
    //vec& operator\(const double s);

  private:
    double* elements;
    int size;

};


vec::vec(int i){
  elements = new double[i];
  size = i;
  for(int k = 0; k < size; ++k ){
    elements[k] = 0.0;
  }
}

vec::~vec(){
  delete[] elements;
}
double& vec::operator[](int i){
  if(i < 0 || i > this->size){
    throw std::runtime_error("Out of bound access");
  }
  return elements[i];
}


const double& vec::operator[](int i)const{

  if(i < 0 || i > this->size){
    throw std::runtime_error("Out of bound access");
  }
  return elements[i];

}
vec vec::operator+(const vec& rhs){
  if(this->get_size() != rhs.get_size()){
      throw std::runtime_error("Incompatible");
  }
  int size = this->get_size();
  vec result(size);
  for(int i = 0; i < size; ++i){
    result[i] = (*this)[i] + rhs[i];
  }

  return result;
}

vec vec::operator-(const vec& rhs){
  if(this->get_size() != rhs.get_size()){
      throw std::runtime_error("Incompatible");
  }
  int size = this->get_size();
  vec result(size);
  for(int i = 0; i < size; ++i){
    result[i] = (*this)[i] - rhs[i];
  }

  return result;

}
vec vec::operator*(const double s){

  int size = this->get_size();
  vec result(size);
  for(int i = 0; i < size; ++i){
    result[i] = (*this)[i] * s;
  }

  return result;   
}
/*vec& vec::operator/(const double s){

  int size = this->get_size();
  vec result(size);
  for(int i = 0; i < size; ++i){
    result[i] = (*this)[i] / s;
  }

  return result;   
}
*/
vec& vec::operator=(const vec& other){
  if(this == &other){
    return *this;
  }

  int size = this->size;
  if(size != other.get_size()){
    throw std::runtime_error("Invalid Dimension in assignment");

  }

  delete[] elements;

  elements = new double[size];

  for(int i = 0; i < size; ++i){
     elements[i] = other[i]; 
  }
  return *this;
}
double dot(const vec& vec_1, const vec& vec_2){
  if(vec_1.get_size() != vec_2.get_size()){
    throw std::runtime_error("Invalid Dimension in dot()");
  }
  double result = 0;
  int size = vec_1.get_size();
  for(int i = 0; i < size; ++i){
    result += vec_1[i] * vec_2[i]; 
  }
  return result;
}

double magnitude(const vec& v){
  int size = v.get_size();
  double result = 0.0;
  for(int i = 0; i < size; ++i){
    result += ( v[i] * v[i]);
  }
  return sqrt(result);
}

int vec::get_size() const{
  return size;
}

void print_vec(const vec& v){
  for(int i =0; i < v.get_size(); ++i){
    printf("%f ", v[i]);    
  
  }
  printf("\n");
}
#endif
