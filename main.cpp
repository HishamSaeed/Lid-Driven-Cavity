#include "helper.hpp"
#include "visual.hpp"
#include "init.hpp"
#include "sor.hpp"
#include "boundary_val.hpp"
#include "uvp.hpp"
#include <cstdio>
#include <iostream>




/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed. Use the predefined matrix<typename> type and give initial values in the constructor.
 * - perform the main loop
 * - at the end: destroy any memory allocated and print some useful statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two-dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop, the following steps are required (for some of the 
 * operations, a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */

int main(int argn, char** args){

  // Acquiring paramters file from user and stop if no file provided
  std::string szFileName;

  // Full path has to be provided in the terminal, not relative
  
  if(argn == 2){
    szFileName = args[1];
    std::cout << "File path read" << szFileName << std::endl;
  }else
  {
    std::cout << "Usage: ./cppfile InputFile (with absolute)" << std::endl;
    return -1;
  }
  

  //Statically allocated no need for dynamic allocations for simulation parameters

  // Solvers Parameters, 
  int itermax;
  double eps;
  double alpha;
  double omg;
  double tau;
  

  // Grid Data
  double xlength;
  double ylength;
  double dx;
  double dy;
  int imax;
  int jmax;

  // Fluid Data, reynolds number, initial conditions for velocity(x and y direction) and pressure
  double Re;
  double UI;
  double VI;
  double PI;

  // Forces Data, gravitational forces
  double GX;  
  double GY;             
  
  // Time Stepping parameters
  double t_end;                
  double dt;             
  double dt_value;

  // Loop parameters
  double t = 0.0;      // time variable for time grid points
  double res = 1.0;    // residual value for the SOR iteration
  int n = 0;           // counter for the variables value at each step
  int it = 0;          // iteration index for the SOR solver

  // Solver required matrices
  matrix<double> F;
  matrix<double> G;
  matrix<double> RS;

  printf("\x1B[96m Reading Simulation parameters from the data file\033[0m\n");      
  // Reading parameters to the static local variables from the parameters file provided by the user
  read_parameters(szFileName,&Re,&UI,&VI,&PI,&GX,&GY,&t_end,&xlength,&ylength,&dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&dt_value);


  // Initialize Grid
  Grid domain = Grid(imax,jmax,1,PI,UI,VI) ;

  // Allocating memory for matrices
  F.resize(domain.imaxb(),std::vector<double>(domain.jmaxb(),0.0));
  G.resize(domain.imaxb(),std::vector<double>(domain.jmaxb(),0.0));
  RS.resize(domain.imaxb(),std::vector<double>(domain.jmaxb(),0.0));
  
  // Setting the boudnary conditions of the grid
   boundaryvalues(imax,jmax,domain);
  
  // 2D Arrays necessary for writing the vtk file
  //Dynamic allocation for the matrices holding the values for grid points
  double** U = new double*[imax+2];
  double** V = new double*[imax+2];
  double** P = new double*[imax+2];

  //Dynamic allocation for the rows of the matrices
  for(int i = 0; i <= imax+1; i++){
    U[i] = new double[jmax+2];
    V[i] = new double[jmax+2];
    P[i] = new double[jmax+2];
  }
  
  printf("\x1B[35m Solving Lid Driven cavity with size %.2fx%.2f, Grid: %dx%d Reynolds: %.2f\033[0m\n",xlength,ylength,imax,jmax,Re);
  printf("\x1B[35m Solver Parameters: omg: %.2f alpha: %.2f dt: %.2f  \033[0m\n",omg,alpha,dt);         

  // Main Loop
  while(t < t_end){

    printf("\x1B[37m t = %.3f n = %d dt = %.5f\033[0m\n",t,n,dt); 

    // Calculating suitable step size
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,domain);

    // Set boundary values for u and v
    boundaryvalues(imax,jmax,domain);
      
    // Compute F n and G n
    calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,domain,F,G);
    // Compute the right-hand side rhs of the pressure equation
    calculate_rs(dt,dx,dy,imax,jmax,F,G,RS);
    
    it = 0;
    res = 1.0;
    // Solve system using SOR
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS,&res);
      it++;
    }

    if(it == itermax)
      printf("\x1B[33m Warning SOR exited due to maximum iteration reached at step: %d\033[0m\n",n);  
      

    calculate_uv(dt,dx,dy,imax,jmax,domain,F,G);
    t = t + dt;
    n = n + 1;
    
  }
  
  printf("\x1B[32m Simulation Finished Succesfully in %d time steps\033[0m\n",n);  
  
  // Getting the velocity and pressure for the vtk file
  domain.velocity(V,velocity_type::V);
  domain.velocity(U,velocity_type::U);
  domain.pressure(P);

  
  // Output values for visualization
  write_vtkFile("Lid_Driven_cavity",n,xlength,ylength,imax,jmax,dx,dy,U,V,P);

  printf("\x1B[32m vtk file written successfully\033[0m\n"); 

  // Free data structures
  free(U);
  free(V);
  free(P);

  return -1;
}
