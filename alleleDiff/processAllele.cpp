//
//  processAllele.cpp
//  alleleDiff
//
//  Created by Gang Peng on 4/16/20.
//  Copyright © 2020 Gang Peng. All rights reserved.
//


#include <iostream>
#include <fstream>

#include "processAllele.hpp"
#include "normal.h"
#include "populationInfo.h"

using namespace std;

bool processAllele(const char* alleleFile, const char* outputFile){
    ifstream fin(alleleFile);
    if(!fin.is_open()){
        cerr << "Cannot ope file "<<alleleFile<<endl;
        return false;
    }
    
    vector<double> maf(population.size(), 0.0);
    int numRace = 0;
    vector<double> mac(population.size(), 0.0);
    vector<double> nchrobs(population.size(), 0.0);
    string currentChr = "";
    string currentSNP = "";
    string currentA1 = "";
    string currentA2 = "";
    
    while (!fin.eof()) {
        string fline;
        string SNP = "";
        getline(fin, fline); // remove header
        getline(fin, fline);
        
        if(fline.size() < 2){
            // process the last record
            break;
        }
        
        vector<string> vsline = split(fline, "\t");
        if(currentSNP != vsline[1]){
            if(currentSNP != ""){
                string outline = currentChr + "\t" + currentSNP + "\t" + currentA1 + "\t" + currentA2 + "\t";
                for(size_t i=0; i<(maf.size()-1); i++){
                    for(size_t j=(i+1); j<maf.size(); j++){
                        outline = outline + to_string(maf[i] - maf[j]) + "\t";
                    }
                }
                
                double macEAS = mac[1] + mac[2] + mac[3] + mac[4] + mac[5];
                double macEUR = mac[6] + mac[7] + mac[8] + mac[9] + mac[10];
                double macAFR = mac[11] + mac[12] + mac[13] + mac[14] + mac[15] + mac[16] + mac[17];
                double macAMR = mac[18] + mac[19] + mac[20] + mac[21];
                double macSAS = mac[22] + mac[23] + mac[24] + mac[25] + mac[26];
                
                double nchrobsEAS = nchrobs[1] + nchrobs[2] + nchrobs[3] + nchrobs[4] + nchrobs[5];
                double nchrobsEUR = nchrobs[6] + nchrobs[7] + nchrobs[8] + nchrobs[9] + nchrobs[10];
                double nchrobsAFR = nchrobs[11] + nchrobs[12] + nchrobs[13] + nchrobs[14] + nchrobs[15] + nchrobs[16] + nchrobs[17];
                double nchrobsAMR = nchrobs[18] + nchrobs[19] + nchrobs[20] + nchrobs[21];
                double nchrobsSAS = nchrobs[22] + nchrobs[23] + nchrobs[24] + nchrobs[25] + nchrobs[26];
                
                double mafEAS = macEAS/nchrobsEAS;
                double mafEUR = macEUR/nchrobsEUR;
                double mafAFR = macAFR/nchrobsAFR;
                double mafAMR = macAMR/nchrobsAMR;
                double mafSAS = macSAS/nchrobsSAS;
                
                
                
            }
        } else {
            
        }
        
    }
    return true;
}
