//
//  complexity.cpp
//  
//
//  Created by Roberto Moran Tovar on 24.05.21.
//
//Template to run a random walk in the antigen space and check for lost of memory

#include "../lib/Immuno_functions.hpp"

#include <stdio.h>
#include <numeric>


int main(int argc, char* argv[]) //argv has 0:L ; 1:L_alphabet ; 2:N_ensemble ; 3:N_epitopes_max ; 4:N_esemble
{
    string Text_files_path = "../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/Complexity/";
    cout<<">Running simulation of Antigen random walk ..."<< endl;
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(r, time(NULL));
    clock_t t1,t2;
    t1=clock();
    //-----------------------------------------------------------------------------
    //Parameters:
    
    int L  = atoi(argv[1]); //length of the sequence
    int L_alphabet = atoi(argv[2]); //length of the alphabet
    int N_epitopes_max = atoi(argv[3]);
    long long int N_ensemble = atoi(argv[4]); // Ensemble size
    long long int T = 150; //Random walk time
    
    //------------Energy Matrix------------------------------------------------------
    vector < vector < double > > MJ;
    MJ.resize(L_alphabet);
    for (int k= 0; k<L_alphabet; k++)
    {
        (MJ[k]).resize(L_alphabet);
    };

    ifstream file("../Input_files/MJ2.txt");
    
    for (unsigned int i = 0; i < L_alphabet; i++) {
        for (unsigned int j = 0; j < L_alphabet; j++) {
            file >> MJ[i][j];
        }
    }

    //------------ Alphabet ----------------------------------------------------------
    //Array with the Alphabet
    vector < string > Alphabet;
    Alphabet.resize(L_alphabet);
    ifstream file2("../Input_files/Alphabet_"+std::to_string(L_alphabet)+".txt");
    cout << "The Alphabet is :";
    for (int k = 0; k < L_alphabet; k++) {

        file2 >> Alphabet[k];
        cout << Alphabet[k] ;
    
    }
    cout << "\n";
    
    vector <vector <int> > Epitopes;
    vector <vector <int> > Abs;
    vector <int> Binding_epi_ab;
    vector <double> State;
    State.resize(T);
    

    double state;
    
    for (int N_epitopes = 1; N_epitopes<=N_epitopes_max; N_epitopes++) {
        
        Epitopes.resize(N_epitopes);
        Abs.resize(N_epitopes);
        Binding_epi_ab.resize(N_epitopes);
        
        for (int n = 0; n<N_epitopes; n++) {
            (Epitopes[n]).resize(L);
            (Abs[n]).resize(L);
        }
        
        cout << N_epitopes << " epitopes..." << endl;
        
        for (int t = 0; t<T; t++) {
            State[t] = 0;
        }
        for (int m = 0; m<N_ensemble; m++) {
            for (int n = 0; n<N_epitopes; n++) {
                // initializing the epitopes with random sequences
                for (int l= 0; l<L; l++)
                {
                    Epitopes[n][l] = randIX(0,L_alphabet-1);
                    //cout << Alphabet[Epitopes[n][l]];
                };
                //cout << endl;
                
                // initializing the abs
                find_complementary(L, L_alphabet, MJ, Epitopes[n], Abs[n]);
                for (int l= 0; l<L; l++)
                {
                    //cout << Alphabet[Abs[n][l]];
                };
                //cout << endl;
                //cout << endl;
                Binding_epi_ab[n] = round(1/(1+exp(-Energy(L, L_alphabet, MJ, Abs[n], Epitopes[n], "MJ", r)-(4.5*L))));
                
            }
            
            state = 1 - std::accumulate(begin(Binding_epi_ab), end(Binding_epi_ab), 1, std::multiplies<double>());
            
            for (int t = 0; t<T; t++) {
                for (int n = 0; n<N_epitopes; n++) {
                    mutate_sequence(L, L_alphabet, Epitopes[n]);
                    Binding_epi_ab[n] = round(1/(1+exp(-Energy(L, L_alphabet, MJ, Abs[n], Epitopes[n], "MJ", r)-(4.5*L))));
                }
                state = 1 - std::accumulate(begin(Binding_epi_ab), end(Binding_epi_ab), 1, std::multiplies<double>());
                State[t] += state;
            }
        }
        ofstream fout (Text_files_path+"state_L-"+std::to_string(L)+"_L_alphabet-"+std::to_string(L_alphabet)+"_N_epitopes-"+ std::to_string(N_epitopes)+".txt");
        
        for (int t = 0; t<T; t++) {
            //cout << State[t]/N_ensemble << "\t";
            fout << State[t]/N_ensemble << "\t";
        }
        
        cout << endl;
        fout.close();
    }
        
    
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;

    return 0;
}
