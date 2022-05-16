#include "uvp.hpp"
#include "helper.hpp"
#include "datastructures.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cfloat>

// ----------------------------------------------------------------------------------------------------
// Determines the value of F and G
void calculate_fg(
        double Re,
        double GX,
        double GY,
        double alpha,
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        Grid& grid,
        matrix<double> &F,
        matrix<double> &G)
{
        // For the sake of readability: declaration/definition of variables needed in calculation of F and G;
        double ui_j,  uip1_j,uim1_j,  ui_jp1,ui_jm1, uim1_jp1; // uip1_j = u[i+1][j], uim1_j = = u[i-1][j]
        double vi_j,  vip1_j,vim1_j,  vi_jp1,vi_jm1, vip1_jm1; // vi_jp1 = v[i][j+1]
        double Dx = 1/dx;
        double Dy = 1/dy;
        double kvis = 1/Re; // kinematic viscosity

        // Extract Velocity matrices from Grid object 
        // Thes matrices can be used to access the velocities instead of calling the cell get function
        matrix<double> umatrix;
        matrix<double> vmatrix;

        grid.velocity(umatrix, velocity_type::U);
        grid.velocity(vmatrix, velocity_type::V);
       
       // Calculation of F --------------------------------------------------------------------------------------
        for (int i=1; i<imax; i++)
        {
                for (int j=1; j<=jmax; j++)
                {
                        // Reading the velocity in the x-direction of the current cell and the sourrounding ones.
                        ui_j    = umatrix[i][j];
                        uip1_j  = umatrix[i+1][j];
                        uim1_j  = umatrix[i-1][j];
                        ui_jp1  = umatrix[i][j+1];
                        ui_jm1  = umatrix[i][j-1];

                        uim1_jp1= umatrix[i-1][j+1];

                        // Reading the velocity in the x-direction of the current cell and the sourrounding ones.
                        vi_j    = vmatrix[i][j];
                        vip1_j  = vmatrix[i+1][j];
                        vim1_j  = vmatrix[i-1][j];
                        vi_jp1  = vmatrix[i][j+1];
                        vi_jm1  = vmatrix[i][j-1];

                        vip1_jm1= vmatrix[i+1][j-1];

                        // Calculation of F[i][j]:
                        F[i][j] = ui_j +
                                  dt*(
                                          // FD Aprroximation of U second order derivative in x and y directions
                                          ( kvis*( (Dx*Dx* (uip1_j-(2*ui_j)+uim1_j) ) + (Dy*Dy* (ui_jp1-(2*ui_j)+ui_jm1) ) ) ) -
                                          // d(u^2)/dx term in Doner scheme:
                                          ( (Dx* ( (pow(0.5*(ui_j+uip1_j),2)) - (pow(0.5*(uim1_j+ui_j),2))) ) + (alpha*Dx* ( ( (0.5*abs(ui_j+uip1_j)) * (0.5*(ui_j-uip1_j)) ) - ( (0.5*abs(uim1_j+ui_j)) * (0.5*(uim1_j-ui_j)) ) )) ) -                                                                       
                                          // d(uv)/dy term in Doner scheme:
                                          ( (Dy* ( ( (0.5*(vi_j+vip1_j)) * (0.5*(ui_j+ui_jp1)) ) - ((0.5*(vi_jm1+vip1_jm1)) * (0.5*(ui_jm1+ui_j))) ) ) + (alpha*Dy* ( ( (0.5*abs(vi_j+vip1_j))*(0.5*(ui_j-ui_jp1)) ) - ( (0.5*abs(vi_jm1+vip1_jm1))*(0.5*(ui_jm1-ui_j)) ) ) ) ) +
                                          // Volume forces
                                          GX
                                     );                      
                        
                }
        }


        // Calculation of G --------------------------------------------------------------------------------------
        for (int i=1; i<=imax; i++)
        {
                for (int j=1; j<jmax; j++)
                {
                        // Reading the velocity in the x-direction of the current cell and the sourrounding ones.

                        ui_j    = umatrix[i][j];
                        uip1_j  = umatrix[i+1][j];
                        uim1_j  = umatrix[i-1][j];
                        ui_jp1  = umatrix[i][j+1];
                        ui_jm1  = umatrix[i][j-1];

                        uim1_jp1= umatrix[i-1][j+1];


                        // Reading the velocity in the x-direction of the current cell and the sourrounding ones.

                        vi_j    = vmatrix[i][j];
                        vip1_j  = vmatrix[i+1][j];
                        vim1_j  = vmatrix[i-1][j];
                        vi_jp1  = vmatrix[i][j+1];
                        vi_jm1  = vmatrix[i][j-1];

                        vip1_jm1= vmatrix[i+1][j-1];

                        // Calculation of G[i][j]:
                        G[i][j] = vi_j +
                                  dt*(
                                          // FD Aprroximation of V second order derivative in x and y directions
                                          ( kvis*( (Dx*Dx* (vip1_j-(2*vi_j)+vim1_j) ) + (Dy*Dy* (vi_jp1-(2*vi_j)+vi_jm1)) ) ) -
                                          // d(uv)/dx term in Doner scheme:
                                          ( (Dx* ( ( (0.5*(ui_j+ui_jp1)) * (0.5*(vi_j+vip1_j)) ) - ((0.5*(uim1_j+uim1_jp1)) * (0.5*(vim1_j+vi_j))) ) ) + (alpha*Dx* ( ( (0.5*abs(ui_j+ui_jp1))*(0.5*(vi_j-vip1_j)) ) - ( (0.5*abs(uim1_j+uim1_jp1))*(0.5*(vim1_j-vi_j)) ) ) ) ) - 
                                          // d(v^2)/dy term in Doner scheme:
                                          ( (Dy* ( (pow(0.5*(vi_j+vi_jp1),2)) - (pow(0.5*(vi_jm1+vi_j),2)) ) ) + (alpha*Dy* ( ( (0.5*abs(vi_j+vi_jp1)) * (0.5*(vi_j-vi_jp1)) ) - ( (0.5*abs(vi_jm1+vi_j)) * (0.5*(vi_jm1-vi_j))) )) ) +
                                          // volume forces
                                          GY
                                     );
                        
                }
        }

        // Updating the Boundary Conditions for F;
        for (int j=1; j<=jmax; j++)
        {
                F[0][j]    = umatrix[0][j];
                F[imax][j] = umatrix[imax][j]; 
        }

        for (int i=1; i<=imax; i++)
        {
                G[i][0]    = vmatrix[i][0];
                G[i][jmax] = vmatrix[i][jmax];
        }

}
// ----------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------
// This operation computes the right hand side of the pressure poisson equation.
void calculate_rs(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &F,
        matrix<double> &G,
        matrix<double> &RS)
{

        for (int i=1; i<=jmax; i++)
        {
                for (int j=1; j<=jmax; j++)
                {
                        RS[i][j] = (1/dt)*( ((F[i][j]-F[i-1][j]) / dx) + ((G[i][j]-G[i][j-1]) / dy));
                }
        }
}

// ----------------------------------------------------------------------------------------------------
// Determines the maximal time step size
void calculate_dt(
        double Re,
        double tau,
        double *dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        Grid &grid)
{
        // Extract Velocity matrices from Grid object 
        // Thes matrices can be used to access the velocities instead of calling the cell get function
        matrix<double> umatrix;
        matrix<double> vmatrix;

        grid.velocity(umatrix, velocity_type::U);
        grid.velocity(vmatrix, velocity_type::V);

        double umax_abs=0, vmax_abs=0;
        double Dx = 1/dx;
        double Dy = 1/dy;

        // Search for the maximum absolute value 
        for (int i =1; i<=imax; i++ )
        {
                for (int j=1; j<=jmax; j++)
                {
                        if (abs(umatrix[i][j]) > umax_abs) { umax_abs = abs(umatrix[i][j]);}
                        if (abs(vmatrix[i][j]) > vmax_abs) { vmax_abs = abs(vmatrix[i][j]);}
                }
        }

        // Checking for the minimum value based on stability conditions
        *dt = tau * std::min( (Re/2)*pow((Dx*Dx)+(Dy*Dy),-1) , std::min( dx/umax_abs , dy/vmax_abs ));
        
}
// ----------------------------------------------------------------------------------------------------
void calculate_uv(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        Grid& grid,
        matrix<double> &F,
        matrix<double> &G)
        
{

        matrix<double> pmatrix;
        grid.pressure(pmatrix);

        // Calculate u
        for (int i =1; i<imax; i++ )
        {
                for (int j=1; j<=jmax; j++)
                {
                        grid.cell(i,j).velocity(velocity_type::U) = F[i][j] - ((dt/dx) * (pmatrix[i+1][j] - pmatrix[i][j]));
                        
                }
        }

        // Calculate v
        for (int i =1; i<=imax; i++ )
        {
                for (int j=1; j<jmax; j++)
                {
                        grid.cell(i,j).velocity(velocity_type::V) = G[i][j] - ((dt/dy) * (pmatrix[i][j+1] - pmatrix[i][j]));
                        
                }
        }


}
// ----------------------------------------------------------------------------------------------------