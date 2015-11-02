#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <lib.h>
//#include <omp.h>
using namespace std;

//Those methods/variables will be used
double random(int &seed);
void initialize(int size, int** grid, bool isRandom, long idum, double& M, double& E);
void readInput(double& tempStart, double& tempMax, double& tempStep, int& size, int& mccycles, bool& isRandomSetup);
void monteCarlo(int** spinArray, int size, long& idum, double* energyDeltas, double& M, double& E);
void output(int size, int cycles, double temp, double* averages);
ofstream ofile;

//main program
int main()
{
    /*
     * PART 1:
     * Initializing everything
     */
    //declare variables
    int size, mccycles;
    long idum;
    double tempStart, tempMax, tempStep, M, E;
    double energyChanges[5], averages[5];
    bool isRandomSetup;
    //Time for different seed every time
    time_t start, finish, seed;
    seed = time(NULL);
    idum = -((long)seed);
    //for output file...
    ofile.open("output.txt");
    ofile << "size" << "\t" << "cycles" << "\t" << "temp" << "\t" << "E" << "\t" << "E²" << "\t" << "M" << "\t" << "M²" << "\t" << "|M|" << endl;
    //read input from screen, see function
    readInput(tempStart, tempMax, tempStep, size, mccycles, isRandomSetup);
    //set up array with spins; 'matrix' is from Morten's lib.cpp
    int** spinArray = (int**) matrix(size, size, sizeof(int));
    //initialize the Grid
    initialize(size, spinArray, isRandomSetup, idum, M, E);
    double temp;

    /*
     * PART 2:
     * performing Monte Carlo
     */
    //loop over desired temperature range
    for (temp=tempStart; temp<=tempMax; temp+=tempStep){
        //reset Energy and magnetization (averages)
        for(int i=0; i<5; i++){
            averages[i]=0;
        }
        //Set array with possible energy changes according to temperature
        for(int i=0; i<5; i++){
            double delEnergy = (4*i)-8;
            energyChanges[i] = exp(((double)-delEnergy)/((double)temp));
        }
        //actual Monte Carlo happens here
        for(int i=0; i<=mccycles; i++){
            monteCarlo(spinArray, size, idum, energyChanges, M, E);
            averages[0]+=E;
            averages[1]+=E*E;
            averages[2]+=M;
            averages[3]+=M*M;
            averages[4]+=abs((int)M);
        }
        //output of data for this temperature
        output(size, mccycles, temp, averages);
    }

    /*
     * PART 3:
     * Cleaning up...
     */
    //free allocated memory
    free_matrix((void **) spinArray);
    //close output
    ofile.close();
    return 0;
}//end of main

//This method takes an integer array and its size and can initialize it in random or ordered way with -1 and 1 (=spins)
void initialize(int size, int** grid, bool isRandom, long idum, double& M, double& E)
{
    int i, j;
    if(isRandom){
        //initialize randomly
        for(i=0; i<size; i++){
            for(j=0; j<size; j++){
                if(ran3(&idum)>0.5){grid[i][j]=1;}else{grid[i][j]=1;}
            }}}else{
        //all spins up!
        for(i=0; i<size; i++){
            for(j=0; j<size; j++){
                grid[i][j]=1;
            }
        }
    }
    M=E=0;
    for(i=0; i<size; i++){
        for(j=0; j<size; j++){
            M+=grid[i][j];
            E-=grid[i][j]*(grid[i][(size+j+1)%size]+grid[(size+i+1)%size][j]);
        }
    }
    return;
}

//This methods reads input about all necessary parameters from screen
void readInput(double& tempStart, double& tempMax, double& tempStep, int& size, int& mccycles, bool& isRandomSetup){
    cout << "Please type in the start temperature!" << endl;
    cin >> tempStart;
    cout << "Please type in the maximum temperature!" << endl;
    cin >> tempMax;
    cout << "Please type in the temperature step length!" << endl;
    cin >> tempStep;
    cout << "Please type in the number of spins per direction!" << endl;
    cin >> size;
    cout << "Please type in the number of Monte Carlo cycles" << endl;
    cin >> mccycles;
    int random;
    cout << "Do you wish a random (1) or ordered (0) setup?" << endl;
    cin >> random;
    if(random==1){isRandomSetup=true;}else{isRandomSetup=false;}
    return;
}

//This is the actual Monte Carlo method as described in the report!
void monteCarlo(int** spinArray, int size, long& idum, double* energyDeltas, double& M, double& E){
    int count;
    for(count=0; count<=(size*size); count++){
        //Pick random position
        int x, y;
        x=(int)ran3(&idum);
        y=(int)ran3(&idum);
        //check energy difference
        double deltaE =(double) 2*spinArray[x][y]*(spinArray[(size+x+1)%size][y]+spinArray[(size+x-1)%size][y]+spinArray[x][(size+y+1)%size]+spinArray[x][(size+y-1)%size]);
        //compare it to random number
        if(ran3(&idum)<=energyDeltas[(int)(deltaE+8)/4]){
            //change spin
            spinArray[x][y]*=-1;
            //update energy and magnetization
            M+=2*spinArray[x][y];
            E+=deltaE;
        }
    }
    return;
}

//This method prints the data of the Monte Carlo method for on setup to previously selected output file
void output(int size, int cycles, double temp, double* averages){
    ofile << size << "\t" << cycles << "\t" << temp << "\t" << averages[0]/cycles << "\t" << averages[1]/cycles << "\t" << averages[2]/cycles << "\t" << averages[3]/cycles << "\t" << averages[4]/cycles << endl;
    return;
}
