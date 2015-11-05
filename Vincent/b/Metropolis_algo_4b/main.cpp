#include <iostream>
#include <math.h>
#include <random>
#include <time.h>
#include <fstream>

#define k_b 1.38064852e-23
#define eV 1.602e-19

using namespace std;

//Performs the Metropolis algorithm
void metropolis_energy(int n, double J, int T, int** s,double* meanM, double* susceptibility, double* meanE);

//Calculates the energy for one specific state
double E(double J, int s12,int s21, int s11, int s22)
{
    return -2.0*J*(double(s12)+double(s21))*(double(s11)+double(s22));
}

//Gives back a random spin value
int spin_generator(int random)
{
    return random%2 == 1 ? 1 : -1;
}

//Generates a random state of the system
void spin_assigner(int** s)
{
    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            s[i][j] = spin_generator(rand());
        }
    }
    return;
}

//Gives the probability for a new state to be accepted
double probability(double dEnergy,int T)
{
    return exp(dEnergy*eV/(k_b*double(T)));
}

int main()
{
    int n;
    double J;
    int Tmin, Tstep, Tmax;
    double meanM = 0;
    double susceptibility = 0;
    double meanE = 0;
    int** s = new int* [2];
    for(int i = 0; i < 2; i++)
    {
        s[i] = new int[2];
    }

    Tmin = 10000;
    Tstep = 100;
    Tmax = 100000;
    n = 5000000;
    J = 1;

    ofstream data("messdaten.txt");

/*
    cout << "The mean energy is: " << meanE << endl
         << "The mean magnetization is: " << meanM << endl
         << "The susceptibility is: " << susceptibility << endl;
*/
    data << "# Temperature, Energy expectationvalue, M expectationvalue, Susceptibility expectationvalue\n";

    for(int T = Tmin; T <= Tmax; T += Tstep)
    {
        metropolis_energy(n, J, T, s, &meanM, &susceptibility, &meanE);

        data << T << "  " << meanE << "    " << meanM << "   " << susceptibility << endl;
    }

    data.close();


    for (int i = 0; i < 2; i++)
    {
        delete[] s[i];
    }
    delete[] s;
    return 0;
}

//This function gives the mean energy value in eV, using the metropolis algorithm
void metropolis_energy(int n, double J, int T, int** s, double* meanM, double* susceptibility, double* meanE)
{
    double Energy, dEnergy, sumEnergies = 0, p, prob;
    double absM, sumabsM = 0, sumM = 0, sumMsquared = 0;
    int x,y, oldspin;

    //initializing the random number generator
    srand(time(NULL));
    default_random_engine generator(rand());
    uniform_real_distribution<double> distribution(0,1);

    //generates the spins
   // spin_assigner(s);

    //Test
    s[0][0] = 1;
    s[0][1] = 1;
    s[1][0] = 1;
    s[1][1] = 1;

    //computes energy to our generated spins
    Energy = E(J, s[0][1],s[1][0],s[0][0],s[1][1]);
    prob = exp(-8*J*eV/(k_b*double(T)));

    cout << "Probability: " << prob << endl << endl;

    //Here the metropolis algorithm starts
    for(int i = 0; i < n ; i++)
    {
        //change to a new state
        x = rand()%2;
        y = rand()%2;
        oldspin = s[x][y];
        s[x][y] *= -1;

        //calculating the difference between the old and new energy
        dEnergy = -4*J*(s[(x+1)%2][y]+s[x][(y+1)%2])*oldspin;

        //saves probability that a step gets accepted
        if(dEnergy < 0)
        {
            p = prob;
        }
        else
        {
            p = 1;
        }

        //decides now, if the step gets accepted
        if(distribution(generator) <= p)
        {
            Energy -= dEnergy;
        }
        else
        {
            s[x][y] = oldspin;
        }

        sumEnergies += Energy;

        absM = abs(s[0][0] + s[0][1] + s[1][0] + s[1][1]);
        sumabsM += absM;
        sumM += s[0][0]+s[0][1]+s[1][0]+s[1][1];
        sumMsquared += absM*absM;

    }

    *meanM = sumabsM/n;
    *susceptibility = (sumMsquared/n-(sumM/n)*(sumM/n))/(k_b*double(T));
    *meanE = sumEnergies/n;

    return;
}

