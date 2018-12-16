#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "random.hpp"

using Random = effolkronium::random_static;
using namespace std;

vector<string> words;
int numFiles = 20;
int docLength = 50;

void crearFiles(){
        for(int i = 0; i < numFiles; i++){
                string path ="./docs/file" + to_string(i) + ".txt";
                ofstream outfile (path);

                for(int j = 0; j < docLength; j++){
                        if(j != 0) outfile << " ";
                        int ran = Random::get() % 50 + 1;
                        outfile << words[ran];
                }
                outfile.close();
        }                
}

int main ()
{          
        cout << "Assegurat de tindre el directopry 'docs' creat" << endl << endl;

        cout << "Introdueix el número de documents que vols crear: ";             
        cin >> numFiles ;
        cout << endl;

        cout << "Introdueix el numero de parauless que vols que tingui cada document: ";
        cin >> docLength;
        cout << endl;

       ifstream input( "paraules.txt" );
       for(string line; getline( input, line ); )
        {
                line.erase(remove(line.begin(), line.end(), '\n'), line.end());
                line.erase(remove(line.begin(), line.end(), '\r'), line.end());
                words.push_back(line);
        }

        crearFiles();
}