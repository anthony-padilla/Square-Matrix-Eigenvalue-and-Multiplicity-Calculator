#ifndef MATRIX_H
#define MATRIX_H
#include<string>
#include<iostream>
#include "config.h"

class matrx{
  
  public:
    //Constructors & Descructor
    matrx(int dim_x, int dim_y);
    matrx();
    ~matrx();
    int  rows() const;
    int columns() const;
    
    double* operator[](const int &i);
    const double* operator[](const int &i)const;
    matrx(const matrx& other);
    matrx& operator=(const matrx& rhs);
    matrx operator*(const matrx& rhs);
    matrx operator*(const double& scalar);
    double** matrix;
  private:
    int dim_x, dim_y;
   

};


matrx::matrx(int x, int y){
  matrix = new double*[x];
  for(int i = 0; i < x; ++i){
    matrix[i] = new double[y];
    for(int j = 0; j < y; ++j){
      matrix[i][j] = 0.0;
    }
  }

  dim_x = x;
  dim_y = y;
}
matrx::~matrx(){
  #ifdef DEBUG
    printf("[DEBUG] Starting Clean\n"); 
  #endif
  for(int i = 0; i < dim_x; ++i){
    delete[] matrix[i];
  } 

  delete[] matrix;
  #ifdef DEBUG
    printf("[DEBUG] cleaned\n");
  #endif
}
int matrx::rows() const{
  return dim_x;
}
int matrx::columns() const{
  return dim_y;
}

double* matrx::operator[](const int &i){
  
  return matrix[i];
}
const double* matrx::operator[](const int &i) const{
  
  return matrix[i];
}

matrx::matrx(const matrx& other):dim_x(other.dim_x), dim_y(other.dim_y){
  matrix = new double*[dim_x];
  for(int i = 0; i < dim_x; ++i){
    matrix[i] = new double[dim_y];
    for(int j = 0; j < dim_y; ++j){
	matrix[i][j] = other.matrix[i][j];
    }
  }
}
matrx& matrx::operator=(const matrx& rhs){
  
  if(this == &rhs){
    return *this;
  }
  for(int i = 0; i < dim_x; ++i){
    delete[] matrix[i];
  }
  delete[] matrix;

  dim_x = rhs.dim_x;
  dim_y = rhs.dim_y;
  
  matrix = new double*[dim_x];

  for(int i = 0; i < dim_x; ++i){
    matrix[i] = new double[dim_y];
    for(int j = 0; j < dim_y; ++j){
	matrix[i][j] = rhs.matrix[i][j];
    }
  }
  return *this;
}

matrx matrx::operator*(const matrx& rhs){
  if (dim_y != rhs.rows()){
    throw std::runtime_error("Cannot Multiply - Dimension Mismatch");
  }

  int rows = dim_x;
  int cols = rhs.columns();
  matrx result(rows, cols);

  for(int i = 0; i < rows; ++i){
    for(int j = 0; j < cols; ++j){
	double sum = 0.0;
	for(int k = 0; k < dim_y;++k){
	  sum+= (*this)[i][k] * rhs[k][j];
	  result[i][j] = sum;
	}
  
    }
  }

  return result;
}
matrx matrx::operator*(const double& scalar){
  
  matrx result(dim_x,dim_y);

  for(int i = 0; i < dim_x; ++i){
    for(int j = 0; j < dim_y;++j){
      result[i][j] = (*this)[i][j] * scalar;
    }
  }

  return result;
}


//## HELPER FUNCTIONS ##



void print_matrix(const matrx &m){
    for(int i = 0; i < m.rows(); i++){
	for(int j = 0; j < m.columns(); j++){
	    double element = m[i][j];
	    if(element >= 0 ){
		    std::cout << element << "  ";
	    }else{
		    std::cout << element << " ";
	    }
	}
	std::cout << std::endl;
    }
    #ifdef DEBUG
      printf("[DEBUG] finished print_matrix()\n");
    #endif
}
void parseArray(std::string inputString,double array[], const int &m){
   
    int i = 0;
    int count = 0;
    std::string num = "";

    while ( count < m){
     
    if(inputString[i] != ' ' && i < inputString.size()){
		num += inputString[i];
	}else {
	
     array[count] = (double) stoi(num);
	 ++count;
	 num = "";
	
    }
	
    ++i;
    }


}

void fill_matrix(matrx &matrix){
  
  const int n = matrix.rows();
  const int m = matrix.columns();

  double* row = new double[m];
  #ifdef DEBUG
    printf("\033[1m fill_array() \033[0m called...\n");
  #endif
  std::string input_string;
  for(int i = 0; i < n; ++i){
    std::getline(std::cin,input_string);
    parseArray(input_string, row, m);
    for(int j = 0; j < m; ++j){
      matrix.matrix[i][j] = row[j];
    }
  }
  delete [] row;
}





#endif
