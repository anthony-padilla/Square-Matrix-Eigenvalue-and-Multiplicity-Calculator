#include<cmath>
#include<iostream>
#include "matrix.h"
#include "linear_algebra.h"
#include "min_pol.h"
#include "config.h"



using namespace std;

int main(){
  
  bool modified = false;

  printf("\033[1m Algebraic and Geometric Multiplicity finder\033[0m \n");
  printf("This program assumes a square n x n matrix\nPlease enter size n: \n");
  int n;
  cin >> n;

  /*MATRIX INPUT*/
  matrx m(n,n);
  printf("[PROGRAM] Dimension n = %d\n", n);
  printf("[PROGRAM] Please Enter Matrix row by row - \n");
  
  cin.ignore(); 
  
  fill_matrix(m);
  printf("Working with matrix - \n");
  print_matrix(m);
  
  /*PART 1 : Isomorphic Matrix*/

  printf("\n[STEP 1] Creating Isomorphic Matrix\n");
  matrx iso_matrix = Isomorphic(m);
  printf("[*] printing Isomorphic() result.\n");
  print_matrix(iso_matrix);

  /*PART 2: Gaussian Elimination on 'iso_matrix' to get minimal polynomial*/
  printf("\n[STEP 2] Reducing Isomorphic Matrix\n");
  matrx reduced_iso_matrix = gaussian_elimination(iso_matrix);
  printf("\n[*] printing gaussian_elimination result...\n");
  print_matrix(reduced_iso_matrix);
  
  /*PART 3: Extract Minimal Polynomial*/
  printf("\n[STEP 3] Extracting Minimal Polynomial\n");
  Min_Pol mp;   
  mp = extract_minimal_polynomial(reduced_iso_matrix);
  printf("[*] completed extract_minimal_polynomial()\n");
  printf("[*] calling print_minimal_polynomial()...\n");

  print_minimal_polynomial(mp);
  
  /*PART 4: Make Companion Matrix*/
  printf("\n[STEP 4] Creating Companion matrix\n");
  matrx comp_matrx = companion_matrix(mp);
  printf("[*] calling print_matrix() for companion matrix...\n");
  print_matrix(comp_matrx);
  
  /*PART 5: Root Finding */
  printf("\n[STEP 5] Root Finding\n");
  printf("[*] Calling find_roots()\n");
  matrx A = find_roots(comp_matrx,modified);
  printf("[*] printing root matrix A\n");
  print_diagonal(A);

  /*PART 6: Determinant*/

  printf("\n[STEP 6] Characterisitc Polynomial and roots\n");
  Min_Pol det_pol;
  printf("[*] Calling fadeev_verrier()\n");
  det_pol = fadeev_verrier(m);
  printf("[*] printing polynomial...\n");
  print_minimal_polynomial(det_pol);
  matrx det_comp_matrx = companion_matrix(det_pol);
  printf("[*] printing companion matrix for char_polynomial..\n");
  print_matrix(det_comp_matrx);
  matrx B = find_roots(det_comp_matrx,modified);
  print_diagonal(B);

  printf("\n\033[1m Results \033[0m\n");
  printf("\033[1m Minimal Polynomial Roots\033[0m\n");
  print_diagonal(A);
  printf("\033[1m Characteristic Polynomial Roots\033[0m\n");
  print_diagonal(B);
  return 0;


}
