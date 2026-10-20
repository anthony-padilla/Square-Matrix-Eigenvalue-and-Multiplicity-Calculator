#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H
#include "matrix.h"
#include "min_pol.h"
#include "vec.h"
matrx Idendity(const int &n){
  
  matrx m(n,n);
  for(int i = 0; i < n; ++i){
      m.matrix[i][i] = 1.0;
  }

  return m;
  
}

/* ##### DESCRIPTION ##### */
/*PARAMETERS: Takes in 2D Matrix object 'matrx', which contains elements, as well as dimension
  OUTPUT: Returns a matrix with
  rows = n^2
  columns = n+ 1
  Each column is the isomorphism of a power of the matrix i.e
  Column 0: A^0
  Colum 1: A^1*/
matrx Isomorphic(matrx m){
  if (m.rows() != m.columns()){
    printf("[ERROR] in Isomprhic() - input is not square\n");
    exit(-1);
  }

  int iso_row_size = m.rows() * m.columns();
  int iso_column_size = m.rows() + 1;
  matrx iso_matrix(iso_row_size, iso_column_size);

  /*Get Identity Easily*/
  matrx identity = Idendity(m.rows());
  for(int i = 0; i < m.rows(); ++i){
   for(int j = 0; j < m.columns(); ++j){
      iso_matrix[m.rows()* i + j ][0] = identity[j][i];
   } 
  }
  matrx power_matrix = m;
  for(int n =1; n < iso_column_size; ++n){
      if(n>1){
	power_matrix = power_matrix * m;
      }
      for(int i = 0; i < m.rows() ; ++i){
	for(int j = 0; j < m.columns(); ++j){
	  iso_matrix[m.rows()*i+j][n] = power_matrix[j][i];
	}
      }

  }

  return iso_matrix;
}

matrx gaussian_elimination( matrx m ){
  #ifdef DEBUG
    printf("[DEBUG] Starting gaussian_elimination()...");
  #endif
  matrx result = m;
    //K is current column but limited by rows 
    
  for(int k = 0 ; k < result.columns(); k++){
      #ifdef DEBUG
	printf("[DEBUG] starting iteration %d \n",k);
      #endif
      int pivot = k;
      double pivot_value = result[pivot][k];
      
      //iterate to find largest pivot for numerical stability
      for(int i = k + 1 ; i < result.rows(); ++i){
	  if(abs(result[i][k]) > pivot_value){
	      
	      pivot_value = result[i][k];
	      pivot = i;

	  }
      }

      if(!result[pivot][k]){
	  continue;
      }

      //Swap rows if we changed pivot
      if(pivot != k){
	  double temp = 0.0;
	  for(int s = 0; s < result.columns(); ++s){
	      temp = result[k][s];
	      result[k][s] = result[pivot][s];
	      result[pivot][s] = temp;
	  }
      }
      //FIXME         
      double normal = result[k][k];
      for(int i = k; i < result.columns(); ++i){
	  if(normal != 0){
	    result[k][i] = result[k][i] / normal;
	  }
      }

      // Zero out column *below* the pivot 
      for(int i = k + 1; i < result.rows(); ++i){
	  double scalar = 0.0;

	  if(result[k][k] > 1e-9){
	   scalar = result[i][k] / result[k][k];
	  }
	  
	  for(int j = k + 1 ; j < result.columns(); ++j){
      
	      result[i][j] -= scalar * result[k][j];
      
	  }
	  
	  result[i][k] = 0; //Done at the end since element in pivot column must be = 0
      }

      //Zero out column *above* the pivot 
      if( k > 0 ){
	for(int i = k -1; i >=0 ; --i){

	  double scalar = 0.0;
	  
	  if(result[k][k] > 1e-9){
	  
	    scalar = result[i][k] / result[k][k];
	  
	  }
	  
	  for(int j = k ; j < result.columns(); ++j){
	  
	    result[i][j] -= scalar * result[k][j];
	 
	  }
	}
      }


    }
  return result;
}

Min_Pol extract_minimal_polynomial(const matrx& m){
  
  bool empty_row = false;
  int i = 0;

  while(!empty_row && i < m.rows()){
      empty_row = true;
      for(int j = 0; j < m.columns(); ++j){
	if(m[i][j]){
	  empty_row =false;
	}
      }
      if(empty_row){
	break;
      }
      ++i;
  }
  if(!empty_row){
    throw std::runtime_error("No Empty Row in extracting minimal polynomial");
  }
  Min_Pol mp(i);
  for(int j = 0; j < i; ++j){
    mp[j] = m[j][i];
  }
  
  return mp;
} 

matrx companion_matrix(const Min_Pol &mp){
  int size = mp.get_size();
  matrx comp_matrx(size,size);
  for(int j = 0; j < size; ++j){
    comp_matrx[0][j] = mp[size-j-1];
  }
  for(int i = 1; i < size; ++i){
    for(int j = 0; j < size; ++j){
      if(j == i-1){
	comp_matrx[i][j] = 1.0;
      }else{
	comp_matrx[i][j] = 0.0;
      }
    }
  }
  return comp_matrx;
}

void gram_schmidt(matrx& A , matrx& Q, matrx& R){
  int n = A.columns();
  for(int j = 0; j < n; ++j){
    vec v_j(n);
    /*Copy NOTE: Could make overload later*/
    for(int k = 0; k < n; ++k){
      v_j[k] = A[k][j];
    }

    for(int i = 0; i < j; ++i){
      vec q_i(n);
      for(int k = 0; k < n; ++k){
	q_i[k] = Q[k][i];
      }
      R[i][j] = dot(q_i,v_j);
      q_i = q_i * R[i][j];
      v_j = v_j - q_i;

    }

    R[j][j] = magnitude(v_j);
    if(R[j][j] < 1e-14){
      R[j][j] = 0.0;
      for(int k =0; k < n; ++k){
	Q[k][j] = 0.0;
      }
      continue;
    }

    v_j = v_j * ( 1.0 / (R[j][j]));
    for(int k = 0; k < n; ++k){
      Q[k][j] = v_j[k];
    }

  }

}

void modified_gram_schmidt(matrx& A, matrx& Q, matrx& R){
  int n = A.columns(); //Could be either, assumes square
  
  for(int k = 0; k < n;++k){
    vec v_k(n);

    for(int i = 0; i < n; ++i){
      v_k[i] = A[i][k];
    }
    
    for(int j=0;j<k;++j){
      vec q_j(n);
      for(int i = 0; i < n; ++i){
        q_j[i] = Q[i][j];
      }

      R[j][k] = dot(q_j,v_k);
      v_k = v_k - q_j * R[j][k];
    }
    R[k][k] = magnitude(v_k);
    v_k = v_k * (1.0/R[k][k]);
    
    for(int i = 0; i < n; ++i){
      Q[i][k] = v_k[i];
    }

  }

}

void qr_iteration(matrx& A, matrx& Q, matrx& R, int max_iterations, bool modified){
  int n = A.columns();
  if(modified){printf("[*] Using Modified GS\n");}else{printf("[*] using Classical GS\n");}
  for(int h = 0; h < max_iterations; ++h){
    for(int i = 0; i < n; ++i){
      for(int j = 0; j < n; ++j){
	  R[i][j] = 0.0;
      }
    }
    if(modified){
      modified_gram_schmidt(A,Q,R);
    }else{
      gram_schmidt(A,Q,R);
    }
    A = R * Q;

    double diag_norm = 0.0;
    for(int r = 1; r < n; ++r){
      for(int c =0; c < r; ++c){
	diag_norm += fabs(A[r][c]);
      }
    }
    if(diag_norm < 1e-8){
      printf("[*] qr_iteration() converged at iteration %d with diag_norm : %.10f\n", h, diag_norm);
      break;
      }
    }
  }


void print_diagonal(const matrx& A){
  for(int i =0; i < A.columns(); ++i){
    printf("%f, ",std::round(A[i][i]));
  }
  printf("\n");
}

/*DESCRIPTION: 
 * PARAMETERS: Takes in matrix A, finds roots
 * OUTPUT: Returns a new matrix A, which contains roots along the diagonal.
 * NOTE: This assumes A is a companion matrix. A is passed-by-copy, so returning A leaves the original untouched (or should)*/
matrx find_roots(matrx A,bool modified){
  int n = A.columns();
  matrx Q(n,n);
  matrx R(n,n);
  
  qr_iteration(A,Q,R,1000000,modified);
  
  return A;
}

/*Assumes square*/
double trace(const matrx& m){
  double sum = 0.0;
  for(int i = 0; i < m.columns(); ++i){
    sum+=m[i][i];
  }
  return sum;
}

Min_Pol fadeev_verrier(matrx& A){
  int n = A.columns();
  Min_Pol mp(n);
  
  matrx M_prev = Idendity(n);

  matrx M_curr(n,n);

  for(int i =1; i <=n; ++i){
    if(i == 1){
      M_curr = A * M_prev;
    }else{
      M_prev = M_curr;

      double t = 0.0;
      t = trace(M_prev);
      
      t = t * 1 / ( i - 1);

      for(int k = 0; k < n; ++k){
	for(int r = 0; r < n; ++r){
	  if(k == r){
	    M_prev[k][r] -= t;
	  }
	}
      }

    M_curr = A * M_prev;

    }
    double c = 0.0;
    c = trace(M_curr);
    c = c * -1.0 * (1.0 / i);
    mp[n-i] = c * -1.0;
  }
  return mp;
}

#endif
