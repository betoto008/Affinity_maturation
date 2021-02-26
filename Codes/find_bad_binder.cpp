//Template to find good binder.
//Input: (all parameters are already set internally!)


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>
#include <algorithm>

using namespace std;
//Library for random number generators
#include "./lib/random.cpp"
//There are two functions extracted from the library
//double randX(min,max): a random number between min and max
//int randIX(min,max): an integer random number between min and max  

//----------------------------------------------------------------------------------

//Function to calculate the energy: Implement the Energy Matrix
double Energy(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen)
{
	double E (0.);

    for(int i=0; i<L ; i++){
    	E = E + MJ[Antigen[i]][sequence[i]];
    }

	return E;
};

//Function to calculate complementary sequence
void find_complementary(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< string > const & Alphabet, vector< int > const & sequence, vector<int> & complementary_sequence)
{
    //cout << "The complementary sequence is: ";
    for(int i=0; i<L ; i++){
        vector < double > v;
        v.resize(L);
        v = MJ[sequence[i]];
        int index = std::min_element(v.begin(), v.end())- v.begin();
        complementary_sequence[i] = index;
        //cout << Alphabet[index];
    }
    //cout << "\n";
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
	cout<<">Running Monte Carlo simulation of the BCRs ...\n"<< endl;
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

	
	//------------ Energy Matrix ------------------------------------------------------
	//MJ Matrix
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

	//------------ Alphabet ----------------------------------------------------------
	//Array with the Alphabet
	vector < string > Alphabet;
	Alphabet.resize(L_alphabet);
	ifstream file2("Alphabet.txt");
	//cout << "The Alphabet is :";

	for (int k = 0; k < L_alphabet; k++) {

	    file2 >> Alphabet[k];
	    //cout << Alphabet[k] ;
	
	}
	//cout << "\n";

	//------------- Initiating Antigen ------------------------------------------------
	//Array with the antigen
	vector < int > Antigen;
	Antigen.resize(L);
	vector < int > Antigen_i;
	Antigen_i.resize(L);

	//---------------Initiating sequence with the Master sequence----------------------
	//Array with the current sequence
	vector < int > master_sequence;
	master_sequence.resize(L);

	//Array with the current sequence
	vector < int > complementary_sequence;
	complementary_sequence.resize(L);
	
	double e = -50;
	double e_new;
	for(int i = 0; i<100000 ; i++){

		for (int k= 0; k<L; k++){
			Antigen_i[k] = randIX(0,L_alphabet-1);
		};

		find_complementary(L, L_alphabet, MJ, Alphabet, Antigen_i, master_sequence);
		e_new = Energy(L, L_alphabet, MJ, master_sequence, Antigen_i);

		if(e_new>e){
			e = e_new;
			Antigen = Antigen_i;
		}
	}
    find_complementary(L, L_alphabet, MJ, Alphabet, Antigen, master_sequence);
    
	cout << "Antigen:";
	for (int k= 0; k<L; k++){
		cout << Alphabet[Antigen[k]];
	};
	cout << "\n";
	cout << "Master Sequence:";
	for (int k= 0; k<L; k++){
		cout << Alphabet[master_sequence[k]];
	};
	cout << "\n";
	cout << "Binding energy:"<< e << "\n";

	/*
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
	fout.close();	
	};

	*/
	
	//------------------------------------------------------------------------------
	cout<< ">Simulation completed…"<< endl;
	t2= clock();
	cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;
	return 0;
}
