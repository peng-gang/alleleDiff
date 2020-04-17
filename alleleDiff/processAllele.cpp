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

bool processAllele(const char* alleleFile, const char* outputFile1, const char* outputFile2){
    ifstream fin(alleleFile);
    if(!fin.is_open()){
        cerr << "Cannot ope file "<<alleleFile<<endl;
        return false;
    }
    
    ofstream fout1(outputFile1);
    ofstream fout2(outputFile2);
    
    fout1 << "CHR\tSNP\tA1\tA2\t";
    for(size_t i = 0; i<(population.size()-1); i++){
        for(size_t j=(i+1); j<population.size(); j++){
            fout1 << population[i] << "-" << population[j] << "\t";
        }
    }
    fout1 <<endl;
    
    fout2 << "CHR\tSNP\tA1\tA2\t";
    for(size_t i = 0; i<(majorPopulation.size()-1); i++){
        for(size_t j=(i+1); j<majorPopulation.size(); j++){
            fout2 << majorPopulation[i] << "-" << majorPopulation[j] << "\t";
        }
    }
    fout2 <<endl;
    
    vector<double> maf(population.size(), 0.0);
    //maf of 5 major population
    vector<double> mafMajor(5, 0.0);
    int numRace = 0;
    vector<double> mac(population.size(), 0.0);
    vector<double> nchrobs(population.size(), 0.0);
    string currentChr = "";
    string currentSNP = "";
    string currentA1 = "";
    string currentA2 = "";
    
    string header;
    getline(fin, header); // remove header
    while (!fin.eof()) {
        string fline;
        getline(fin, fline);
        
        if(fline.size() < 2){
            // process the last record
            if(numRace != population.size()){
                cout << "Race information is missing for " + currentSNP <<endl;
            } else {
                string outline = currentChr + "\t" + currentSNP + "\t" + currentA1 + "\t" + currentA2 + "\t";
                string outline2 = outline;
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
                
                mafMajor[0] = mafEAS;
                mafMajor[1] = mafEUR;
                mafMajor[2] = mafAFR;
                mafMajor[3] = mafAMR;
                mafMajor[4] = mafSAS;
                
                for(size_t i=0; i<(mafMajor.size()-1); i++){
                    for(size_t j=(i+1); j<mafMajor.size(); j++){
                        outline2 = outline2 + to_string(mafMajor[i] - mafMajor[j]) + "\t";
                    }
                }
                
                fout1 << outline << endl;
                fout2 << outline2 << endl;
            }
            break;
        }
        
        vector<string> vsline = split(fline, " ", true);
        if(currentSNP != vsline[1]){
            if(currentSNP != ""){
                if(numRace != population.size()){
                    cout << "Race information is missing for " + currentSNP <<endl;
                    currentChr = vsline[0];
                    currentSNP = vsline[1];
                    currentA1 = vsline[3];
                    currentA2 = vsline[4];
                    continue;
                }
                
                string outline = currentChr + "\t" + currentSNP + "\t" + currentA1 + "\t" + currentA2 + "\t";
                string outline2 = outline;
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
                
                mafMajor[0] = mafEAS;
                mafMajor[1] = mafEUR;
                mafMajor[2] = mafAFR;
                mafMajor[3] = mafAMR;
                mafMajor[4] = mafSAS;
                
                for(size_t i=0; i<(mafMajor.size()-1); i++){
                    for(size_t j=(i+1); j<mafMajor.size(); j++){
                        outline2 = outline2 + to_string(mafMajor[i] - mafMajor[j]) + "\t";
                    }
                }
                
                fout1 << outline << endl;
                fout2 << outline2 << endl;
                numRace = 0;
            }
            currentChr = vsline[0];
            currentSNP = vsline[1];
            currentA1 = vsline[3];
            currentA2 = vsline[4];
        }
        
        int idx = -1;
        for(size_t i=0; i<population.size(); i++){
            if(population[i] == vsline[2]){
                idx = (int)i;
                break;
            }
        }
        
        if(idx < 0){
            cout << "Cannot find population " + vsline[2] + " in database" <<endl;
            fin.close();
            fout1.close();
            fout2.close();
            return false;
        }
        
        maf[idx] = stod(vsline[5]);
        mac[idx] = stod(vsline[6]);
        nchrobs[idx] = stod(vsline[7]);
        numRace++;
        
    }
    
    fin.close();
    fout1.close();
    fout2.close();
    
    return true;
}
