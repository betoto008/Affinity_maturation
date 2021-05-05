//
//  Immuno_functions.hpp
//  
//  Created by Roberto Moran Tovar on 27.02.21.
//

#ifndef Immuno_functions_h
#define Immuno_functions_h


#endif /* Immuno_functions_h */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>

//Library for random number generators
#include "./random.cpp"
//There are two functions extracted from the library
//double randX(min,max): a random number between min and max
//int randIX(min,max): an integer random number between min and max

const long double N_A = 6.02214076E23;

using namespace std;

// ---------------- CLASSES ----------------

class bcell {
public:
    vector < int > seq; //vector with the positions of the aa sequence in the alphabet
    bcell();
    bcell(int const & L, int const & L_alphabet, vector< int > & seq);
    double e; //energy with respect to the current epitope.
    double cs;
    bool plasma;
    bool GC;
    bool active;
};
bcell::bcell(){

}
bcell::bcell(int const & L, int const & L_alphabet, vector< int > & seq){
    this->seq = seq;
    cs = 1.0;
    plasma = 0;
    GC = 0;
    active = 0;
    
}

// ---------------- FUNCTION ---------------

//Function to calculate the energy: Implement the Energy Matrix
double Energy(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen)
{
    double E (0.);

    for(int i=0; i<L ; i++){
        E = E + MJ[Antigen[i]][sequence[i]];
    }
    return E;
};

//Function to calculate the energy difference due to a mutation
inline double delt( int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen, int const & pos, int const & aa)
{
    double deltaE (0.);
    deltaE = MJ[Antigen[pos]][aa] - MJ[Antigen[pos]][sequence[pos]];
    return deltaE;
};

//Function to change from aminoacids to positions
void aa_to_positions( int const & L, int const & L_alphabet, vector< string > & Alphabet,  vector< int > & sequence_pos,  string  sequence_aa)
{
    for(int i=0; i<L ;i++){
        
        for(int j=0; j<L_alphabet ;j++){
            if(sequence_aa[i] == Alphabet[j][0]){
                sequence_pos[i] = j;
            }
        }
    }
};

//Function to calculate complementary sequence
void find_complementary(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector<int> & complementary_sequence)
{
    //cout << "The complementary sequence is: ";
    for(int i=0; i<L ; i++){
        vector < double > v;
        v.resize(L);
        v = MJ[sequence[i]];
        int index = std::min_element(v.begin(), v.end()) - v.begin();
        complementary_sequence[i] = index;
    }
    //cout << "\n";
};

// Function to return energies < E_max

void high_energy_filter(){
}

// Function to generate the initial amount of sequences
void generate_Bcells(int N, int L, int L_alphabet, vector<bcell> & Bcells){
    
    //---------Array with the current Sequence-------------------------------------------
    vector < int > Sequence;
    Sequence.resize(L);
    
    for(int n =0 ; n<N ; n++){
        
        //Initiating Sequence with random sequence------------------------------------
        for (int k= 0; k<L; k++)
        {
            Sequence[k] = randIX(0,L_alphabet-1);
        };
        
        // Create a bcell and add it to the vector
        bcell bcell_i(L, L_alphabet, Sequence);
        Bcells[n]  = bcell_i;
    }
}

// Function that selects the antigen-specific naive Bcells from all the sequences
void choose_naive_Bcells(int N, int L, int L_alphabet, vector< vector<double> > const & MJ, vector< int > const & Antigen, vector<bcell> & Bcells, vector<bcell*> & Naive, int & n_naive){
    
    vector <int> MS;
    MS.resize(L);
    find_complementary(L, L_alphabet, MJ, Antigen, MS);
    double e_MS = Energy(L, L_alphabet, MJ, MS, Antigen);
    double e;
    for(int n = 0 ; n<N ; n++){
        e = Energy(L, L_alphabet, MJ, Bcells[n].seq, Antigen);
        if(e<e_MS+28){
            Bcells[n].e = e;
            Naive.push_back( &Bcells[n]);
            n_naive++;
        }
    }
}

//Function that mutates any sequences in one random position
void mutate_sequence(int L, int L_alphabet, vector< int > & sequence){
    int pos = randIX(0,L-1);
    int aa = randIX(0,L_alphabet-1);
    if(sequence[pos] != aa){
        sequence[pos]=aa;
    }
}



