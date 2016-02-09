#include <smoothField/Core>
#include <smoothField/Tree>
#include <smoothField/Interpolate>
#include <iostream>

int main(){

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> coord(81,2);
  coord << -2, -2,
    -1.5, -2,
    -1, -2,
    -0.5, -2,
    0, -2,
    0.5, -2,
    1, -2,
    1.5, -2,
    2, -2,
    -2, -1.5,
    -1.5, -1.5,
    -1, -1.5,
    -0.5, -1.5,
    0, -1.5,
    0.5, -1.5,
    1, -1.5,
    1.5, -1.5,
    2, -1.5,
    -2, -1,
    -1.5, -1,
    -1, -1,
    -0.5, -1,
    0, -1,
    0.5, -1,
    1, -1,
    1.5, -1,
    2, -1,
    -2, -0.5,
    -1.5, -0.5,
    -1, -0.5,
    -0.5, -0.5,
    0, -0.5,
    0.5, -0.5,
    1, -0.5,
    1.5, -0.5,
    2, -0.5,
    -2, 0,
    -1.5, 0,
    -1, 0,
    -0.5, 0,
    0, 0,
    0.5, 0,
    1, 0,
    1.5, 0,
    2, 0,
    -2, 0.5,
    -1.5, 0.5,
    -1, 0.5,
    -0.5, 0.5,
    0, 0.5,
    0.5, 0.5,
    1, 0.5,
    1.5, 0.5,
    2, 0.5,
    -2, 1,
    -1.5, 1,
    -1, 1,
    -0.5, 1,
    0, 1,
    0.5, 1,
    1, 1,
    1.5, 1,
    2, 1,
    -2, 1.5,
    -1.5, 1.5,
    -1, 1.5,
    -0.5, 1.5,
    0, 1.5,
    0.5, 1.5,
    1, 1.5,
    1.5, 1.5,
    2, 1.5,
    -2, 2,
    -1.5, 2,
    -1, 2,
    -0.5, 2,
    0, 2,
    0.5, 2,
    1, 2,
    1.5, 2,
    2, 2; 

  Eigen::Matrix<double,Eigen::Dynamic,1> obs(81);
  obs << -16, -11.375, -9, -8.125, -8, -7.875, -7, -4.625, 0, -11.375, -6.75, -4.375, -3.5, -3.375, -3.25, -2.375, 0, 4.625, -9, -4.375, -2, -1.125, -1, -0.875, 0, 2.375, 7, -8.125, -3.5, -1.125, -0.25, -0.125, 0, 0.875, 3.25, 7.875, -8, -3.375, -1, -0.125, 0, 0.125, 1, 3.375, 8, -7.875, -3.25, -0.875, 0, 0.125, 0.25, 1.125, 3.5, 8.125, -7, -2.375, 0, 0.875, 1, 1.125, 2, 4.375, 9, -4.625, 0, 2.375, 3.25, 3.375, 3.5, 4.375, 6.75, 11.375, 0, 4.625, 7, 7.875, 8, 8.125, 9, 11.375, 16;
  smoothField<double> test(coord,obs);

  smoothTree<double> ttrae(1.5,&test);

  Eigen::Matrix<double,Eigen::Dynamic,1> minCoords(2);
  minCoords << -0.5,-0.5;
  Eigen::Matrix<double,Eigen::Dynamic,1> maxCoords(2);
  maxCoords << 0.5,0.5;

  bicubicInterpolation<double> bci(&test,minCoords,maxCoords);

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> nc(121,2);
  nc << -0.5, -0.5,
    -0.4, -0.5,
    -0.3, -0.5,
    -0.2, -0.5,
    -0.1, -0.5,
    0, -0.5,
    0.1, -0.5,
    0.2, -0.5,
    0.3, -0.5,
    0.4, -0.5,
    0.5, -0.5,
    -0.5, -0.4,
    -0.4, -0.4,
    -0.3, -0.4,
    -0.2, -0.4,
    -0.1, -0.4,
    0, -0.4,
    0.1, -0.4,
    0.2, -0.4,
    0.3, -0.4,
    0.4, -0.4,
    0.5, -0.4,
    -0.5, -0.3,
    -0.4, -0.3,
    -0.3, -0.3,
    -0.2, -0.3,
    -0.1, -0.3,
    0, -0.3,
    0.1, -0.3,
    0.2, -0.3,
    0.3, -0.3,
    0.4, -0.3,
    0.5, -0.3,
    -0.5, -0.2,
    -0.4, -0.2,
    -0.3, -0.2,
    -0.2, -0.2,
    -0.1, -0.2,
    0, -0.2,
    0.1, -0.2,
    0.2, -0.2,
    0.3, -0.2,
    0.4, -0.2,
    0.5, -0.2,
    -0.5, -0.1,
    -0.4, -0.1,
    -0.3, -0.1,
    -0.2, -0.1,
    -0.1, -0.1,
    0, -0.1,
    0.1, -0.1,
    0.2, -0.1,
    0.3, -0.1,
    0.4, -0.1,
    0.5, -0.1,
    -0.5, 0,
    -0.4, 0,
    -0.3, 0,
    -0.2, 0,
    -0.1, 0,
    0, 0,
    0.1, 0,
    0.2, 0,
    0.3, 0,
    0.4, 0,
    0.5, 0,
    -0.5, 0.1,
    -0.4, 0.1,
    -0.3, 0.1,
    -0.2, 0.1,
    -0.1, 0.1,
    0, 0.1,
    0.1, 0.1,
    0.2, 0.1,
    0.3, 0.1,
    0.4, 0.1,
    0.5, 0.1,
    -0.5, 0.2,
    -0.4, 0.2,
    -0.3, 0.2,
    -0.2, 0.2,
    -0.1, 0.2,
    0, 0.2,
    0.1, 0.2,
    0.2, 0.2,
    0.3, 0.2,
    0.4, 0.2,
    0.5, 0.2,
    -0.5, 0.3,
    -0.4, 0.3,
    -0.3, 0.3,
    -0.2, 0.3,
    -0.1, 0.3,
    0, 0.3,
    0.1, 0.3,
    0.2, 0.3,
    0.3, 0.3,
    0.4, 0.3,
    0.5, 0.3,
    -0.5, 0.4,
    -0.4, 0.4,
    -0.3, 0.4,
    -0.2, 0.4,
    -0.1, 0.4,
    0, 0.4,
    0.1, 0.4,
    0.2, 0.4,
    0.3, 0.4,
    0.4, 0.4,
    0.5, 0.4,
    -0.5, 0.5,
    -0.4, 0.5,
    -0.3, 0.5,
    -0.2, 0.5,
    -0.1, 0.5,
    0, 0.5,
    0.1, 0.5,
    0.2, 0.5,
    0.3, 0.5,
    0.4, 0.5,
    0.5, 0.5; 


  std::cout << "Test af LPR:\n";
  std::cout << "newval1 <- c(";
  for(int l = 0; l < nc.rows(); ++l){
    std::cout << test(nc.row(l))(0) << ((l < nc.rows()-1) ? ", " : ") \n\n");
  }

  std::cout << "Test af Tree:\n";
  std::cout << "newval2 <- c(";
  for(int l = 0; l < nc.rows(); ++l){
    std::cout << ttrae(nc.row(l))(0) << ((l < nc.rows()-1) ? ", " : ") \n\n");
  }

  // NOTE: Result from Tree should NOT be the same as for bicubic and cubic (which should be the same) because the points used for interpolation are not the same!
  
  std::cout << "Test af biCubic interpolation:\n";
  std::cout << "newval3 <- c(";
  for(int l = 0; l < nc.rows(); ++l){
    std::cout << bci(nc.row(l))(0) << ((l < nc.rows()-1) ? ", " : ") \n\n");
  }

std::cout << "Test af cubic interpolation:\n";
  cubicInterpolation<double> testci(&test,minCoords,maxCoords);
  std::cout << "newval4 <- c(";
  for(int l = 0; l < nc.rows(); ++l){
    std::cout << testci(nc.row(l))(0) << ((l < nc.rows()-1) ? ", " : ") \n\n");
  }


  
  return 0;
}
