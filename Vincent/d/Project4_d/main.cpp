#include <iostream>
#include <math.h>
#include <random>
#include <time.h>
#include <fstream>

#define k_b 1.38064852e-23
#define eV 1.602e-19

using namespace std;

//Performs the Metropolis algorithm
void metropolis_energy(int n, double J, int T, int** s,double* meanM, double* susceptibility, double* meanE, ofstream& data);

//Calculates the energy for one specific state
double E(double J, int** s)
{
    double sum[3] = {0,0,0};
    for(int i = 0; i < 20; i++)
    {
        for(int j = 0; j < 20; j++)
        {
            sum[0] += s[i][j]*(s[i][(j+1)%20]+s[(i+1)%20][j]);
        }
        sum[1] += s[i][19]*s[i][0];
        sum[2] += s[0][i]*s[19][i];
    }
    return -J*(sum[0]+sum[1]+sum[2]);
}

//Gives back a random spin value
int spin_generator(int random)
{
    return random%2 == 1 ? 1 : -1;
}

//Generates a random state of the system
void spin_assigner(int** s)
{
    for(int i = 0; i < 20; i++)
    {
        for(int j = 0; j < 20; j++)
        {
            s[i][j] = spin_generator(rand());
        }
    }
    return;
}

int main()
{
    int n;
    double J;
    int Tmin, Tstep, Tmax, T;
    double meanM = 0;
    double susceptibility = 0;
    double meanE = 0;
    int** s = new int* [20];
    for(int i = 0; i < 20; i++)
    {
        s[i] = new int[20];
    }

    T = 30000;
    Tmin = 10000;
    Tstep = 100;
    Tmax = 100000;

    n = 100000000;
    J = 1;
    ofstream data("messdaten_d.txt");

    data << "# State; Frequency of a state\n";

    metropolis_energy(n, J, T, s, &meanM, &susceptibility, &meanE, data);

/*
    for(int T = Tmin; T <= Tmax; T += Tstep)
    {
        data << T << "  ";

        metropolis_energy(n, J, T, s, &meanM, &susceptibility, &meanE, data);

        data << meanE << "    " << meanM << "   " << susceptibility << endl;
    }
*/
    data.close();


    for (int i = 0; i < 20; i++)
    {
        delete[] s[i];
    }
    delete[] s;
    return 0;
}

//This function gives the mean energy value in eV, using the metropolis algorithm
void metropolis_energy(int n, double J, int T, int** s, double* meanM, double* susceptibility, double* meanE, ofstream& data)
{
    double Energy, dEnergy, sumEnergies = 0, p, prob[5];
    double M = 0, absM, sumabsM = 0, sumM = 0, sumMsquared = 0;
    int x, y, alpha, oldspin;
    int Energy_state_amount[200], state;

    for(int i = 0;i < 200; i++)
    {
        Energy_state_amount[i] = 0;
    }

    //initializing the random number generator
    srand(time(NULL));
    default_random_engine generator(rand());
    uniform_real_distribution<double> distribution(0,1);

    //generates the spins
    //spin_assigner(s);
    for(int i = 0; i < 20; i++)
    {
        for(int j = 0; j < 20; j++)
        {
            s[i][j] = 1;
        }
    }
    //Every spin up gives the lowest energy. Hence we are in state number 0
    state = 0;

    //calculates magnetization
    for(int i = 0; i < 20 ; i++)
    {
        for(int j = 0; j < 20; j++)
        {
            M += s[i][j];
        }
    }

    //computes energy to our generated spins
    Energy = E(J, s);

    prob[2] = exp(-4*J*eV/(k_b*double(T)));
    prob[4] = exp(-8*J*eV/(k_b*double(T)));
    cout << "Probabilities: " << prob[4] << "     " << prob[2] << endl << endl;

    //Here the metropolis algorithm starts
    for(int i = 0; i < n ; i++)
    {
        //change to a new state
        x = rand()%20;
        y = rand()%20;
        oldspin = s[x][y];
        s[x][y] *= -1;

        alpha = s[(x+1)%20][y]+s[(x+19)%20][y]+s[x][(y+1)%20]+s[x][(y+19)%20];

        //calculating the difference between the old and new energy
        dEnergy = -2*J*alpha*oldspin;

        //saves probability that a step gets accepted
        if(dEnergy < 0)
        {
            p = prob[alpha];
        }
        else
        {
            p = 1;
        }

        //decides now, if the step gets accepted
        if(distribution(generator) <= p)
        {
            Energy -= dEnergy;
            state += (alpha*oldspin)/2;
        }
        else
        {
            s[x][y] = oldspin;
        }

        sumEnergies += Energy;
        M += s[x][y] - oldspin;
        absM = abs(M);
        sumabsM += absM;
        sumM += M;
        sumMsquared += absM*absM;

        Energy_state_amount[state]++;
    }

    *meanM = sumabsM/n;
    *susceptibility = (sumMsquared/n-(sumM/n)*(sumM/n))/(k_b*double(T));
    *meanE = sumEnergies/n;

    for(int i = 0; i < 200; i++)
    {
        data << i << "  " << Energy_state_amount[i] << endl;
    }

    return;
}
