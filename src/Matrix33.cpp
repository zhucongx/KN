#include "Matrix33.h"

// #define ARMA_ALLOW_FAKE_GCC
// #define ARMA_DONT_USE_WRAPPER
// Uncomment this line if there is a compilation error
#include "armadillo"

Matrix33<double> InverseMatrix33(const Matrix33<double> &input) {
  arma::mat mat_input = {{input.row1.x, input.row1.y, input.row1.z},
                         {input.row2.x, input.row2.y, input.row2.z},
                         {input.row3.x, input.row3.y, input.row3.z}};
  arma::mat inverse_matrix = arma::inv(mat_input);
  return {{inverse_matrix(0, 0), inverse_matrix(0, 1), inverse_matrix(0, 2)},
          {inverse_matrix(1, 0), inverse_matrix(1, 1), inverse_matrix(1, 2)},
          {inverse_matrix(2, 0), inverse_matrix(2, 1), inverse_matrix(2, 2)}};
  // double det = (input.row1.x * input.row2.y * input.row3.z
  //     - input.row1.x * input.row2.z * input.row3.y
  //     - input.row1.y * input.row2.x * input.row3.z
  //     + input.row1.y * input.row2.z * input.row3.x
  //     + input.row1.z * input.row2.x * input.row3.y
  //     - input.row1.z * input.row2.y * input.row3.x);
  // return {{(input.row2.y * input.row3.z - input.row2.z * input.row3.y) / det,
  //          (input.row1.z * input.row3.y - input.row1.y * input.row3.z) / det,
  //          (input.row1.y * input.row2.z - input.row1.z * input.row2.y) / det},
  //         {(input.row2.z * input.row3.x - input.row2.x * input.row3.z) / det,
  //          (input.row1.x * input.row3.z - input.row1.z * input.row3.x) / det,
  //          (input.row1.z * input.row2.x - input.row1.x * input.row2.z) / det},
  //         {(input.row2.x * input.row3.y - input.row2.y * input.row3.x) / det,
  //          (input.row1.y * input.row3.x - input.row1.x * input.row3.y) / det,
  //          (input.row1.x * input.row2.y - input.row1.y * input.row2.x) / det}};

}