//Template to run a Monte Carlo simulation for the BCRs.
//Input: (all parameters are already set internally!)

//Output: output.txt (columns: 1 temperature , 2 energy per spin, 3 heat capacity, 4 magnetisation per spin, 5 susceptibility  )

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>

using namespace std;
//Library for random number generators
#include "./lib/random.cpp"
//There are two functions extracted from the library
//double randX(min,max): a random number between min and max
//int randIX(min,max): an integer random number between min and max  

//----------------------------------------------------------------------------------

//Function to calculate the energy: Implement the Energy Matrix
double energy(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen)
{
	double E (0.);

    for(int i=0; i<L ; i++){
    	E = E + MJ[Antigen[i]-1][sequence[i]-1];
    }

	return E;
};

//Function to calculate the energy difference due to a mutation
//Do NOT use energy(seq 1)-energy(seq 2) as this is computatioLally costly: the energy difference depends on mutation position only
inline double delt( int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen, int const & pos, int const & aa)
{
	double deltaE (0.);
	deltaE = MJ[Antigen[pos]-1][aa-1] - MJ[Antigen[pos]-1][sequence[pos]-1];
	return deltaE;
};

//----------------------------------------------------------------------------------
int main(int argc, char* argv[])
{	
	string Text_files_path = "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/";
	cout<<">Running Monte Carlo simulation of the BCRs ..."<< endl;
	clock_t t1,t2;
    t1=clock();
	//-----------------------------------------------------------------------------
	//Parameters: (they are fine as they are)
	int L (9); //length of the sequence
	int L_alphabet (20);
	int nT (1); //Number of temperature points
	double T1 (.1) ; double T2 (2);
	long long int n0 (0*L), N02 (2E6*L), d0 (10*L); //Number of steps: initial prelude, total, distance between sampling points
	int N0[1] = {1E8};
	int array1[9] = {8, 5, 6, 17, 14, 18, 9, 8, 17};
	int array2[9] = {3, 3, 3, 3, 3, 7, 4, 3, 3};
	//Array with the antigen
	vector < int > Antigen;
	Antigen.resize(L);
	//Array with the current sequence
	vector < int > sequence;
	sequence.resize(L);
	//------------Energy Matrix------------------------------------------------------
	vector < vector < double > > MJ;
	MJ.resize(L_alphabet);
	for (int k= 0; k<L_alphabet; k++)
	{
		(MJ[k]).resize(L_alphabet);
	};

	ifstream file("MJ2.txt");

	for (unsigned int i = 0; i < L_alphabet; i++) {
	    for (unsigned int j = 0; j < L_alphabet; j++) {
	        file >> MJ[i][j];
	    }
	}
	//--------------------------------------------------------------------------------
	//Initiating Antigen -------------------------------------------------------------
	for (int k= 0; k<L; k++)
	{
		Antigen[k]= array1[k];
	};

	for (int kT= 0; kT<nT; kT++)
	{	

		//Initiating sequence with the Master sequence------------------------------------
		for (int k= 0; k<L; k++)
		{
			sequence[k]= array2[k];
			//sequence[k] = randIX(1,L_alphabet);
		};
		//--------------------------------------------------------------------------------

		// Set the temperature
		//double T (T2-(T2-T1)*double(kT)/double(nT-1));
		double T = T2;

		//Output file
		ofstream fout (Text_files_path+"output_highT_N-"+ std::to_string(N0[kT])+".txt");

		cout<< ">T= "<< T<< endl;
		
		double E; //Energy
		E= energy(L,L_alphabet,MJ,sequence,Antigen);
		fout<< E<< "\t"<<0.0 << "\t"<< 0 <<endl;
        
        //Starting the trajectory:
		int countData (0); //Number of data point sampled in the trajectory
		for (long long int k= 0; k < (N0[kT]*L); k++)
		{	
			//FOR YOU TO FILL-IN:
			//Pick up a position and an aminoacid and calculate the energy difference if it were mutated
  			int pos = randIX(0,L-1);
  			int aa = randIX(1,L_alphabet);
  			
  			double deltaE = delt(L, L_alphabet,MJ, sequence, Antigen, pos, aa);

			//Decide whether to actually flip the spin or not: (Metropolis' algorithm)
			if (deltaE<0){
				sequence[pos] = aa;
			}
			else{
				double rand = randX(0,1);
				if(rand < exp((-1*deltaE)/T)){
					sequence[pos]=aa;
				}
			};
			//
			
			//Calculate the observables starting after n0 steps: sample data points every d0 steps (for Eq.(2))
			if (k>=n0)
			{
				if ( (k%d0) == 0)
				{
					//Increase the number of data sum
					countData++;
					E= energy(L,L_alphabet,MJ,sequence,Antigen);
					fout<< E<< "\t"<< deltaE << "\t"<< aa <<endl;				
				};
			};
		};
		/*
        //Divide the sums by the number of data points to get the averages
		avrE/=double(countData); avrEE/=double(countData); avrM/=double(countData); avrMM/=double(countData);
		
		//FOR YOU TO FILL IN:
		double heatcapacity= (avrEE - avrE*avrE)/(T*T) ;
        double susceptibility= (avrMM - avrM*avrM)/(T);
		

		//Write to an output file:
		//fout<< T<< "\t"<< avrE/double(N)<< "\t"<< heatcapacity/double(N)<< "\t"<< avrM/double(N)<< "\t"<< susceptibility/double(N) <<endl;
		*/
		fout.close();
	};

	
	//------------------------------------------------------------------------------
	cout<< ">Simulation completed…"<< endl;
	t2= clock();
	cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;
	return 0;
}
