#include "sor.hpp"
#include <cmath>

void sor(
        double omg,
        double dx,
        double dy,
        int    imax,
        int    jmax,
        Grid& grid,
        matrix<double> &RS,
        double *res
) {
    int i,j;
    double rloc;
    double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
    matrix<double> P;
    grid.pressure(P);

    /* SOR iteration */
    for(i = 1; i <= imax; i++) {
        for(j = 1; j<=jmax; j++) {
            P.at(i).at(j) = (1.0-omg)*P.at(i).at(j)
                            + coeff*(( P.at(i+1).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)+P.at(i).at(j-1))/(dy*dy) - RS.at(i).at(j));
        }
    }

    /* compute the residual */
    rloc = 0;
    for(i = 1; i <= imax; i++) {
        for(j = 1; j <= jmax; j++) {
            rloc += ( (P.at(i+1).at(j)-2.0*P.at(i).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)-2.0*P.at(i).at(j)+P.at(i).at(j-1))/(dy*dy) - RS.at(i).at(j))*
                    ( (P.at(i+1).at(j)-2.0*P.at(i).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)-2.0*P.at(i).at(j)+P.at(i).at(j-1))/(dy*dy) - RS.at(i).at(j));
        }
    }
    rloc = rloc/(imax*jmax);
    rloc = sqrt(rloc);
    /* set residual */
    *res = rloc;


    /* set boundary values */
    for(i = 1; i <= imax; i++) {
        P.at(i)[0] = P.at(i)[1];
        P.at(i).at(jmax+1) = P.at(i).at(jmax);
    }
    for(j = 1; j <= jmax; j++) {
        P[0].at(j) = P[1].at(j);
        P.at(imax+1).at(j) = P.at(imax).at(j);
    }

    grid.set_pressure(P);

}

