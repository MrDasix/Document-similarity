#include <iostream>
#include <algorithm>
#include <ctime>
#include <map>
#include <vector>
#include <string>
#include <list>
#include <set>
#include <fstream>
#include <functional>
#include <climits>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem::v1;
using namespace std;

const bool DEBUG = 0;
const int shingleSize = 5 ;

map<int,vector<string> > docsInfo;
map<int, set<string> > docsAsShingleSets_String;
map<int, set<int> > docsAsShingleSets;
set<int> docNames;
int numShingles, numHashes;
vector<int> tp; 
vector<vector<int> > signatures;

vector<string> split(const string& s, char delimiter)
{
	vector<string> tokens;
	string token;
	istringstream tokenStream(s);
	while (getline(tokenStream, token, delimiter))
	{
                if(token.size() > 0 &&  token != " "){
                        token.erase(remove(token.begin(), token.end(), ' '), token.end());
			tokens.push_back(token);     
		} 
	}
	return tokens;
}

void llegirDocuments(){
	string path = "./docs";
	int id = 0;
	docNames.clear();
	for (const auto & entry : fs::directory_iterator(path)){
		docNames.insert(id);

		ifstream t(entry.path());
		string str;

		t.seekg(0, t.end);   
		str.reserve(t.tellg());
		t.seekg(0, t.beg);

		str.assign((istreambuf_iterator<char>(t)),istreambuf_iterator<char>()); //read document as string
		str.erase(remove(str.begin(), str.end(), '\n'), str.end());
		str.erase(remove(str.begin(), str.end(), '\r'), str.end());

		docsInfo[id] = split(str , ' ');		

                if(DEBUG){
                        cout << "Document " <<  id << " has all this words: " << endl;
                        for(int j = 0; j < docsInfo[id].size(); j++){
                                cout << docsInfo[id][j] <<  endl;
                        }
                        cout << endl << endl;
                }
		id++;
	}
}

void convertirShingles(){
	for(int i= 0; i < docNames.size(); i++){
		for(int k =0 ; k < docsInfo[i].size(); k++){
			if(k+shingleSize <= docsInfo[i].size()){			
				string shingle = "";
				for(int j=0; j < shingleSize; j++){
					if(j != 0 ) shingle += " ";
					shingle += docsInfo[i][k+j];
				}
				docsAsShingleSets_String[i].insert(shingle);//string shingle

                                hash<std::string> hasher;
                                unsigned long hashed = hasher(shingle);
                                int abbreviated_hash = hashed & INT_MAX;
                                docsAsShingleSets[i].insert(abbreviated_hash);//int hashed shingle
			}
		}	

                if(DEBUG){
                        cout << "Document " <<  i << " has all this shingles: " << endl;
                        for(auto f : docsAsShingleSets[i]) {
                                cout << f << endl;
                        } 
                        cout << endl << endl;
                }
	}
}

///////////////////      Triangle Matrix

// Guarda els valors de similitud. Per guardar entre parelles, només
// la meitat dels elements de tota la matriu. Usant un Triangle Matrix només
// usem menys de la meitat de la memoria per tota la matriu i protegeix d'accedir
// a cel·les vuides o invalides de la matriu.
int getTriangleIndex(int i, int j, map<int, set<int> >& docsAsShingleSets) {
    if (i == j) {
        cout << "No es pot accedir a la matriu triangle amb i == j" << endl;
        exit(0);
    }
    if (j < i) {
        // Swap dels valors
        int tmp = i;
        i = j;
        j = tmp;
    }
    return int(i * (docsAsShingleSets.size() - (i + 1) / 2.0) + j - i) - 1;
}

// mirar si funciona
struct comp {
    bool operator() (const int& a, const int& b) const {
        return a > b;
    }
};

struct comp2 {
    bool operator() (const int& a, const int& b) const {
        return a > b;
    }
};

void inputJaccardSimilarity(int &docid, int& veins, int numDocs) {
    bool b = false;
    while (!b) {
        cout << "\nIntrodueix la id del document que t'interessa: (Dintre del rang [1-" << numDocs << "])" << endl;
        cin >> docid;
        cout << "input " << docid << endl;
        if (docid <= 0 || docid > numDocs) cout << "\nLa id introduida está fora de rang..." << endl;
        else b = true;
    }
    b = false;
    while (!b) {
        cout << "\nIntrodueix el número de veins propers/similars que vols trobar:" << endl;
        cin >> veins;
        if (veins <= 0) cout << "\nHa de ser un número positiu...";
        else b = true;
    }
}

void JaccardSimilarity(int docid, int veins, map<int, set<int> >& docsAsShingleSets, set<int>& docNames, vector<int>& tp) {
    // Numero d'elements que necessita la matriu Triangle
    int numElems = (docsAsShingleSets.size() * docsAsShingleSets.size() - 1) / 2;

    // JSim els valors reals de Jaccard Similarity
    vector<double> JSim (numElems);

    ///////////////////      Jaccard Similarities

    cout << endl << endl << "Calculant Jaccard Similarities del Shingles..." << endl;

    clock_t t0 = clock();

    map<int, set<int> >::iterator first = docsAsShingleSets.begin();

    int i = docid;

    // com un loading progress cada 100 documents
    //if (i % 100 == 0) cout << "     (" << i << "/" << docAsShingleSets.size() << ")" << endl;

    map<int, set<int> >::iterator iti = docsAsShingleSets.find(i);
    set<int> s1 = docsAsShingleSets[iti->first];
    map<int, int, comp> veinsDelDocI;           // ordenats pel major percentatge

    map<int, set<int> >::iterator itj = docsAsShingleSets.begin();
    while (itj != docsAsShingleSets.end()) {
        if (itj->first != iti->first) {
            // Agafar el conjunt de shingles del document j
            set<int> s2 = docsAsShingleSets[itj->first];
            // Calcula i guarda el Jaccard similarity
            set<int> interseccio;
            set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), /*back_*/inserter(interseccio, interseccio.begin()));

            set<int> unio;
            set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), /*back_*/inserter(unio, unio.begin()));

            JSim[getTriangleIndex(iti->first,itj->first,docsAsShingleSets)] = interseccio.size() / double(unio.size());
            double percsimilar = JSim[getTriangleIndex(iti->first,itj->first,docsAsShingleSets)] * 100;
            if (percsimilar > 0) {
                std::fixed;
                cout.precision(2);
                cout << "\t" << iti->first << "\t --> " << itj->first << "\t" << percsimilar << "%" << endl;
                veinsDelDocI[itj->first] = percsimilar;
            }
        }
        ++itj;
    }

    cout << endl << "Comparant Shingles..." << endl;
    cout << endl << "Els top " << veins << " documents més similars al document " << iti->first << " son: " << endl;

    map<int,int>::iterator iterator= veinsDelDocI.begin();
    while (iterator != veinsDelDocI.end()) {
        tp.push_back(iterator->first);
        cout << endl << "Shingles del Document " << iterator->first << " amb Jaccard Similarity " << iterator->second << "%" << endl;
        iterator++;
    }

    double elapsed_secs = double (clock()-t0) / CLOCKS_PER_SEC;

    cout << endl << "Ha tardat " << elapsed_secs << " segons en calcular totes les Jaccard Similarities de Shingles" << endl;
    //cout << "Tot i que tarda més que amb k-shingles o minhash o lsh, aquests percentatges son els reals" << endl;

    // Esborrar JSim aqui si causa problemes de ram
}

void displayAllSignaturesISimilaritat(vector<vector<int> >& signatures, int docid, int veins, map<int, set<int> >& docsAsShingleSets, vector<int>& tp, int numHashes) {
    cout << endl << "Numero de signatures: " << signatures.size(); // list signatures size

    int numElems = (docsAsShingleSets.size() * docsAsShingleSets.size() - 1) / 2;
    // estJSim els valors estimats de Jaccard Similarity comparant amb MinHash signatures
    vector<double> estJSim (numElems);

    int tpsig = 0;  // true positives
    int fpsig = 0;  // false positives

    clock_t t0 = clock();

    double threshold = 0;
    cout << endl <<"Jaccard Similarity entre signatures" << endl;
    cout << endl <<"Els valors mostrats son els valors estimats de Jaccard Similarity" << endl;

    int i = docid;
    vector<int> signature1 = signatures[i];

    map<int, int, comp2> veinsDelDocI;  // veins ordenats per major percentatge de similitud

    map<int, set<int> >::iterator itj = docsAsShingleSets.begin();

    while(itj != docsAsShingleSets.end()) {
        if (i != itj->first) {
            vector<int> signature2 = signatures[itj->first];
            // contem quantes posicions del minhash son iguals
            int count = 0;
            for (int k = 0; k < numHashes; ++k) {
                if (signature1[k] == signature2[k]) ++count;
            }

            estJSim[getTriangleIndex(i, itj->first, docsAsShingleSets)] = (count / float(numHashes));

            if (estJSim[getTriangleIndex(i,itj->first,docsAsShingleSets)] > 0) {
                set<int> s1 (signature1.begin(), signature1.end());
                set<int> s2 (signature2.begin(), signature2.end());

                set<int> interseccio;
                set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(interseccio, interseccio.begin()));
                set<int> unio;
                set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(unio, unio.begin()));
                double p = (interseccio.size() / float(unio.size()));

                if (double(p) > threshold) {
                    double percsimilar = estJSim[getTriangleIndex(i,itj->first,docsAsShingleSets)] * 100;
                    veinsDelDocI[itj->first] = p * 100;
                }
            }
        }
    }

    vector<int> sigpos;
    cout << endl << "Comparant Signatures" << endl;
    cout << endl <<"Els top " << veins << " més similars al document " << docid << " son:" << endl;

    map<int, int, comp2>::iterator it = veinsDelDocI.begin();
    while (it != veinsDelDocI.end()) {
        cout << endl <<"Signatures del Document " << it->first << " amb Jaccard Similarity " << it->second << "%" << endl;
        sigpos.push_back(it->first);
        ++it;
    }

    set<int> s(tp.begin(), tp.end());
    set<int> inters;
    set_intersection(tp.begin(), tp.end(), sigpos.begin(), sigpos.end(), inserter(inters, inters.begin()));
    fpsig = veins - inters.size();
    tpsig = veins - fpsig;
    double elapsed = clock() - t0;
    cout << endl << tpsig << "/" << veins << " True Positives i " << fpsig << "/" << veins << " False Positives Generats al comparar signatures." << endl;
    cout << endl << "Ha tardat " << elapsed << " secons en calcular Jaccard Similarity amb Signatures." << endl;
}

bool esPrimer(int numero) {
    if(numero < 2)
        return false;
    if(numero == 2)
        return true;
    if(numero % 2 == 0)
        return false;
    for(int i = 3; i <= sqrt(numero); i += 2)
    {
        if(numero % i == 0)
            return false;
    }
    return true;
}

int getPrimerPrimer() {
    int i;
    for(i=1; not esPrimer(numShingles+i); ++i);
    return numShingles+i;
}

bool conte(vector<int> vec, int n){
    for(int i=0; i<vec.size(); ++i){
        if(vec[i] == n) return true;
    }
    return false;
}

vector<int> obteCoeficients(int n) {
    srand(time(0));
    vector<int> coeficients;
    while(coeficients.size() < n){
        int r = rand()%numShingles;
        if(not conte(coeficients, r)) coeficients.push_back(r);
    }
    return coeficients;
}

int main()
{
	llegirDocuments();
	convertirShingles();

        int primerPrimer = getPrimerPrimer();
        vector<int> coeficients1 = obteCoeficients(numHashes), coeficients2 = obteCoeficients(numHashes);
        vector<vector<int>> signatures;
        for(int docId: docNames){
                set<int> shingles = docsAsShingleSets[docId];
                vector<int> signatura;
                for(int i=0; i<numHashes; ++i){
                        int minHashCode = primerPrimer + 1;
                        for(int shingle: shingles){
                                int hashCode = (coeficients1[i] * shingle + coeficients2[i]) % primerPrimer;
                                if(hashCode < minHashCode) minHashCode = hashCode;
                        }
                        signatura.push_back(minHashCode);
                }
                signatures.push_back(signatura);
        }

        int docid;
        int veins;

        inputJaccardSimilarity(docid, veins, docNames.size());
        JaccardSimilarity(docid, veins, docsAsShingleSets, docNames, tp);
        displayAllSignaturesISimilaritat(signatures, docid, veins, docsAsShingleSets, tp, numHashes);
}