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

using namespace std;

//Function to calculate the energy: Implement the Energy Matrix
double energy(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen)
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
void find_complementary(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< string > const & Alphabet, vector< int > const & sequence, vector<int> & complementary_sequence)
{
    //cout << "The complementary sequence is: ";
    for(int i=0; i<L ; i++){
        vector < double > v;
        v.resize(L);
        v = MJ[sequence[i]];
        int index = std::min_element(v.begin(), v.end()) - v.begin();
        complementary_sequence[i] = index;
        //cout << Alphabet[index];
    }
    //cout << "\n";
};
