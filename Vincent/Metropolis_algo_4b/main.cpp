#include <iostream>
#include <math.h>
#include <random>
#include <time.h>

#define k_b 1.38064852e-23

using namespace std;

double metropolis_energy(int n, double J, double T, int** s,double* meanM);

double E(double J, int s12,int s21, int s11, int s22)
{
    return -2.0*J*(double(s12)+double(s21))*(double(s11)+double(s22));
}

int spin_generator(int random)
{
    return random%2 == 1 ? 1 : -1;
}

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

double probability(double Energy, double nEnergy,double T)
{
    return exp((Energy - nEnergy)/(k_b*T));
}

int main()
{
    int n;
    double J, T;
    double meanM;
    int** s = new int* [2];
    for(int i = 0; i < 2; i++)
    {
        s[i] = new int[2];
    }

    T = 10;
    n = 1000;
    J = 1.602e-19;

/*
    srand(time(NULL));
    default_random_engine generator(rand());
    uniform_real_distribution<double> distribution(0,1);
*/

    cout << "The mean energy is: " << metropolis_energy(n, J, T, s, &meanM) << " eV" << endl
         << "The mean magnetization is: " << meanM << endl;





    for (int i = 0; i < 2; i++)
    {
        delete[] s[i];
    }
    delete[] s;
    return 0;
}

//This function gives the mean energy value in eV, using the metropolis algorithm
double metropolis_energy(int n, double J, double T, int** s, double* meanM)
{
    double Energy, nEnergy, sumEnergies = 0, sumM = 0, p;

    //initializing the random number generator
    srand(time(NULL));
    default_random_engine generator(rand());
    uniform_real_distribution<double> distribution(0,1);

    //generates the spins
    spin_assigner(s);
    //computes energy to our generated spins
    Energy = E(J, s[0][1],s[1][0],s[0][0],s[1][1]);

    //Here the metropolis algorithm starts
    for(int i = 0; i < n ; i++)
    {
        spin_assigner(s);
        nEnergy = E(J, s[0][1],s[1][0],s[0][0],s[1][1]);

        //saves probability that a step gets accepted
        if(Energy < nEnergy)
        {
            p = probability(Energy, nEnergy, T);
        }
        else
        {
            p = 1;
        }

        //decides now, if the step gets accepted
        if(distribution(generator) <= p)
        {
            Energy = nEnergy;
        }

        sumEnergies += Energy;
        sumM += abs(s[0][0])+abs(s[0][1])+abs(s[1][0])+abs(s[1][1]);

    }

    *meanM = sumM/n;
    return sumEnergies/(n*1.602176e-19);
}

