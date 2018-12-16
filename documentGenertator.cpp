#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h> 
#include <time.h>   
#include <unistd.h>
#include <stdio.h>
#include <limits.h>
#include <algorithm>
#include <string>

using namespace std;

vector<string> words;
const int numFiles = 20;
const int docLength = 200;

void crearFiles(){
        for(int i = 0; i < numFiles; i++){
                string path ="./docs/file" + to_string(i) + ".txt";
                ofstream outfile (path);

                for(int j = 0; j < docLength; j++){
                        if(j != 0) outfile << " ";
                        int ran = rand() % 50 + 1;
                        outfile << words[ran];
                }
                outfile.close();
        }                
}

int main ()
{                       
       srand (time(NULL));
       ifstream input( "paraules.txt" );
       for(string line; getline( input, line ); )
        {
                line.erase(remove(line.begin(), line.end(), '\n'), line.end());
                line.erase(remove(line.begin(), line.end(), '\r'), line.end());
                words.push_back(line);
        }

        crearFiles();
}