//*********************************************************
// Jianping Xiao
// AOSS, University of Michigan, Ann Arbor
// March-23-2013
// Project: Multigrid Solvers
// Course: MATH 671 (Fast Algorithm)
// Project Duration: ~ 3 weeks
//********************************************************* 

# include <iostream>
# include <fstream>
# include <cmath>
# include <new>
# include <iomanip>
# include "constant.h"
# include "multigrid_header.h"
using namespace std;

int main(int  argc, char* argv[]){

    int          N   = 64;
    double tolfactor = 0.10;
    ofstream myfile;
    myfile.open ("data.txt");

    int Nxy = N+1; 
    cout<<"grid sizes:"<<endl;            
    cout<<"N:  "<<N <<endl;
    cout<<"Nxy: "<<Nxy<<endl;
    
    double xystart   = 0.0; 
    double xyend     = 1.0;
    double dxy       = (xyend-xystart)/(Nxy-1);
    double *xy       = new double[Nxy];
    double *vxy      = new double[Nxy*Nxy];
    double *vxynew   = new double[Nxy*Nxy];
    double *vxyexact = new double[Nxy*Nxy];
    double *fv       = new double[Nxy*Nxy];
    double *vr       = new double[Nxy*Nxy];
    double *fv2      = new double[(N/2+1)*(N/2+1)];

    int    i, j, k, counter, ixy, ixyup, ixydown, ixyleft, ixyright;  
    int    nlevel, nfactor;                   
    double L2error = 1.0, L2residual = 1.0, residual = 1.0; 
    double lambda  = 1.0;
    double fcoef   = 1.0/(dxy*dxy);     //   lambda*lambda/(dxy*dxy);    //     
  
    for(i=0;i<Nxy;i++){
        xy[i] = xystart + i*dxy;
    }
 
    // initialization;
    initialsolver( Nxy, xy, vxy, vxyexact, fv);


    counter = 0;
    while(L2residual>0.000000000005){

         counter = counter + 1;
         vcycle( N, dxy, vxy, vxynew, fv, (int) log2(N) - 3);    // V_cycle;
         L2error = 0.0;
         L2residual = 0.0;
           for(i=1;i<(Nxy-1);i++){
             for(j=1;j<(Nxy-1);j++){
                ixy          = j*Nxy+i;
                ixyup        =  ixy + Nxy;  
                ixydown      =  ixy - Nxy;   
                ixyleft      =  ixy - 1;     
                ixyright     =  ixy + 1; 
                L2error      = L2error + (vxy[ixy]-vxyexact[ixy])*(vxy[ixy]-vxyexact[ixy]) ;                                // for jacobsolver;
                residual     = -fcoef*( vxy[ixyup]+vxy[ixydown]+vxy[ixyleft]+vxy[ixyright] - 4.0*vxy[ixy] ) - fv[ixy];      // for jacobsolver;
                L2residual   =  L2residual + residual*residual;
             }
          }  
         L2residual = sqrt( L2residual/((Nxy-1)*(Nxy-1)) );
         L2error = sqrt( L2error/((Nxy-1)*(Nxy-1)) );
         cout<<counter<<setprecision(15)<<"    L2error:  "<<L2error<<"   L2residual:  "<<L2residual<<endl;
         myfile<<counter<<"     "<<L2error<<"     "<<L2residual<<endl;

    }


}

