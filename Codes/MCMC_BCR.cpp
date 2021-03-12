//
//  MCMC_BCR.cpp
//  
//  Created by Roberto Moran Tovar on 27.02.21.
//
//Template to run a Monte Carlo simulation for the BCRs.
//Input: (all parameters are already set internally!)

#include "./lib/Immuno_functions.hpp"

//Library for random number generators
#include "./lib/random.cpp"
//There are two functions extracted from the library
//double randX(min,max): a random number between min and max
//int randIX(min,max): an integer random number between min and max

//----------------------------------------------------------------------------------
using namespace std;

//----------------------------------------------------------------------------------
int main(int argc, char* argv[]) //argv has 1:L , 2:N , 3:T
{
    string Text_files_path = "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/MCMC/";
    cout<<">Running Monte Carlo simulation of the BCRs ..."<< endl;
    clock_t t1,t2;
    t1=clock();
    //-----------------------------------------------------------------------------
    //Parameters:
    int L  = atoi(argv[1]); //length of the sequence
    int L_alphabet (20); //length of the alphabet
    int NT (1); //Number of temperature runs
    std::string T2_s (argv[3]);
    double T1 (.1) ; double T2  = std::stod(T2_s); //range of temperatures
    long long int n0 (0*L), d0 (100*L); //Number of steps: initial prelude, distance between sampling points
    long N0[1] = {atoi(argv[2])}; // Array with MCMC steps

    //------------Energy Matrix------------------------------------------------------
    vector < vector < double > > MJ;
    MJ.resize(L_alphabet);
    for (int k= 0; k<L_alphabet; k++)
    {
        (MJ[k]).resize(L_alphabet);
    };

    ifstream file("MJ2.txt");

    //------------ Alphabet ----------------------------------------------------------
    //Array with the Alphabet
    vector < string > Alphabet;
    Alphabet.resize(L_alphabet);
    ifstream file2("Alphabet.txt");
    cout << "The Alphabet is :";
    for (int k = 0; k < L_alphabet; k++) {

        file2 >> Alphabet[k];
        cout << Alphabet[k] ;
    
    }
    cout << "\n";
    for (unsigned int i = 0; i < L_alphabet; i++) {
        for (unsigned int j = 0; j < L_alphabet; j++) {
            file >> MJ[i][j];
        }
    }
    //------------- Antigen -------------------------------------------------------------
    //Array with the antigen
    string Antigen_aa;
    cout << "Insert the Aminoacid sequence of the antigen:\n";
    getline(cin, Antigen_aa);
    
    vector < int > Antigen;
    Antigen.resize(L);
    aa_to_positions(L, L_alphabet, Alphabet, Antigen, Antigen_aa);

    //-----------------------------------------------------------------------------------
    
    //---------Array with the Master Sequence--------------------------------------------
    vector < int > Master_Sequence;
    Master_Sequence.resize(L);
    find_complementary(L, L_alphabet, MJ, Antigen, Master_Sequence);
    double E0 = Energy(L, L_alphabet, MJ, Master_Sequence, Antigen);
    
    //---------Array with the current Sequence-------------------------------------------
    vector < int > Sequence;
    Sequence.resize(L);
    Sequence = Master_Sequence;
    for (int k= 0; k<L; k++)
    {
        cout << Alphabet[Master_Sequence[k]];
    }
    cout << "\n";
    
    //initialize current energy of the MCMC
    double E;
    E= Energy(L,L_alphabet,MJ,Sequence,Antigen);
    cout << E <<"\n";
    
    // For-loop over different temperatures
    for (int kT = 0; kT<NT; kT++)
    {
        //Initiating Sequence with random sequence------------------------------------
        //for (int k= 0; k<L; k++)
        //{
            //Sequence[k] = randIX(1,L_alphabet);
        //};
        //--------------------------------------------------------------------------------

        // Set the temperature
        //double T (T2-(T2-T1)*double(kT)/double(nT-1));
        double T = T2;

        //Output file
        ofstream fout (Text_files_path+"output_L-"+std::to_string(L)+"_T-"+std::to_string(T)+"_N-"+ std::to_string(N0[kT])+"_Antigen-"+Antigen_aa+".txt");

        cout<< ">T= "<< T<< endl;
        
        fout<< E <<endl;
        
        // Starting the MCMCM trajectory:
        int countData (0); //Number of data point sampled in the trajectory
        std::cout << "N= " << N0[kT] << std::endl;
        for (long long int k= 0; k < (N0[kT]*L); k++)
        {
            //Pick up a position and an aminoacid and calculate the energy difference if it were mutated
              int pos = randIX(0,L-1);
              int aa = randIX(0,L_alphabet-1);
              
              double deltaE = delt(L, L_alphabet, MJ, Sequence, Antigen, pos, aa);

            //Decide whether to actually mutate the sequence: (Metropolis' algorithm)
            if (deltaE<0){
                Sequence[pos] = aa;
            }
            else{
                double rand = randX(0,1);
                if(rand < exp((-1*deltaE)/T)){
                    Sequence[pos]=aa;
                }
            };
            //Calculate and print the observables starting after n0 steps: sample data points every d0 steps
            if (k>=n0)
            {
                if ( (k%d0) == 0)
                {
                    E= Energy(L,L_alphabet,MJ,Sequence,Antigen);
                    //if (E<(E0+12)) {//if the energy is lower than E0+12
                        // print the energy value
                        fout<< E << endl;
                    //}
                }
            }
        }
        
        fout.close();
    };

    //------------------------------------------------------------------------------
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;
    return 0;
}

