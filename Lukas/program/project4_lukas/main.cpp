#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <lib.h>
using namespace std;
//the following fuctions should be used in the full program:
double random(int &seed);
void initialize(int latticesize, int ** Spinarray,double &M, double &E,double temp);
void readinfromscrean(int &latticesize, double &Tempstart, double &Tempstep, double &Tempmax);
void metropolis(int latticesize, long &idum , int **Spinarray,double &E, double &M, double *w );
void outputfunction(int latticesize, int maxcycle, double temp, double* average);

/*      Begin of definition of different functions
 *
 *
 *
*/
void initialize(int latticesize, int **Spinarray, double &M, double &E,double temp){
    for(int i=0;i<latticesize;i++){
        for(int j;j<latticesize;j++){
            Spinarray[i][j]=1;
        }
    }
    M=0;
    E=0;
    for(int i=0;i<latticesize;i++){
        for(int j;j<latticesize;j++){
            M+=Spinarray[i][j];
            E-=Spinarray[i][j]*(Spinarray[(latticesize+1+i)%latticesize][j]+Spinarray[i][(latticesize+j+1)%latticesize]);
        }
    }
}

void readinfromscrean(int &latticesize, double &Tempstart, double &Tempstep, double &Tempmax){
    cout << "give me a latticesize please";
    cin >> latticesize;
    cout << "give me a start- Temperature";
    cin >> Tempstart;
    cout << "give me a Temperaturestep";
    cin >> Tempstep;
    cout << "give me a Maximumtemperature";
    cin >> Tempmax;
}

void metropolis(int latticesize, long &idum , int **Spinarray,double &E, double &M, double *w ){
    for(int i=0;i<latticesize;i++){
        for(int j=0;j<latticesize;j++){
            //we now flip one spin at a time, therefore, we have to have random numbers
            int ix= (int) (ran3(&idum)*(double)latticesize);
            int iy= (int) (ran3(&idum)*(double)latticesize);
            int DeltaE;
            DeltaE=2*Spinarray[ix][iy]*(Spinarray[(ix+1)%latticesize][iy]+Spinarray[(ix-1)%latticesize][iy]+Spinarray[ix][(iy+1)%latticesize]+Spinarray[ix][(iy-1)%latticesize]);
            //Do we want to accept this move??
            if(ran3(&idum)<=w[DeltaE+8]){
                Spinarray[ix][iy]*=-1;
                //calculate new Energy and Magnetisation
                M+=(double) 2*Spinarray[ix][iy];
                E+=(double) DeltaE;
            }
        }
    }
}

void outputfunction(int latticesize, int maxcycle, double temp, double* average){

}


/*      Begin of the main function
 *
 *
 *
*/
int main(){
    time_t seed;
    long idum;
    seed=time(NULL);
    int latticesize,maxcycle=100;
    double Magnetisation,Energy,averagevalues[5];
    double Tempstart,Tempstep,Tempmax,E,M,w[17];
    readinfromscrean(latticesize,Tempstart,Tempstep,Tempmax);
    int **grid= new int *[latticesize];
    for(int i=0;i<latticesize;i++){
        grid[i]= new int [latticesize];
    }

    idum=-1;

    for(double temp=Tempstart;temp<Tempmax;temp+=Tempstep){
        E=M=0;
        //set up prob for Energys
        for(int i=0;i<18;i++){
            w[i]=0;
        }
        for(int i=0;i<18;i=i+4){
            w[i]=exp((-8+i)/temp);
        }
        //setting averagevalues to 0;
        for(int i=0;i<=5;i++){
            averagevalues[i]=0.0;
        }
        initialize(latticesize,grid,Magnetisation,Energy,temp);
        //now performing Monte Crlo method
        metropolis(latticesize,idum,grid,E,M,w);
        //compute averages:
        for(int cycles=1;cycles<maxcycle;cycles++){
            averagevalues[0]+=E;
            averagevalues[1]+=E*E;
            averagevalues[2]+=M;
            averagevalues[3]+=M*M;
            averagevalues[4]+=fabs(M);
        }
    outputfunction(latticesize, maxcycle, temp, averagevalues);
    }
}

