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

using namespace std;

std::vector<std::string> population {
    //EAS
    "CHB",
    "JPT",
    "CHS",
    "CDX",
    "KHV",
    //EUR
    "CEU",
    "TSI",
    "FIN",
    "GBR",
    "IBS",
    //AFR
    "YRI",
    "LWK",
    "GWD",
    "MSL",
    "ESN",
    "ASW",
    "ACB",
    //AMR
    "MXL",
    "PUR",
    "CLM",
    "PEL",
    //SAS
    "GIH",
    "PJL",
    "BEB",
    "STU",
    "ITU"
};

std::vector<std::string> EAS {
    "CHB",
    "JPT",
    "CHS",
    "CDX",
    "KHV"
};

std::vector<std::string> EUR {
    "CEU",
    "TSI",
    "FIN",
    "GBR",
    "IBS"
};

std::vector<std::string> AFR {
    "YRI",
    "LWK",
    "GWD",
    "MSL",
    "ESN",
    "ASW",
    "ACB"
};

std::vector<std::string> AMR {
    "MXL",
    "PUR",
    "CLM",
    "PEL"
};

std::vector<std::string> SAS {
    "GIH",
    "PJL",
    "BEB",
    "STU",
    "ITU"
};


std::vector<std::string> majorPopulation {
    "EAS",
    "EUR",
    "AFR",
    "AMR",
    "SAS"
};

bool processAllele(const char* alleleFile, const char* outputFile1, const char* outputFile2,
                   double mafCutoff){
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
    size_t numRace = 0;
    vector<double> mac(population.size(), 0.0);
    vector<double> nchrobs(population.size(), 0.0);
    string currentChr = "";
    string currentSNP = "";
    string currentA1 = "";
    string currentA2 = "";
    string currentTag = "";
    
    string header;
    getline(fin, header); // remove header
    size_t numSNP = 0;
    while (!fin.eof()) {
        string fline;
        getline(fin, fline);
        
        if(fline.size() < 2){
            // process the last record
            if(numRace != population.size()){
                cout << "Race information is missing for " + currentSNP <<endl;
            } else {
                double macEAS = mac[0] + mac[1] + mac[2] + mac[3] + mac[4];
                double macEUR = mac[5] + mac[6] + mac[7] + mac[8] + mac[9];
                double macAFR = mac[10] + mac[11] + mac[12] + mac[13] + mac[14] + mac[15] + mac[16];
                double macAMR = mac[17] + mac[18] + mac[19] + mac[20];
                double macSAS = mac[21] + mac[22] + mac[23] + mac[24] + mac[25];
                
                double nchrobsEAS = nchrobs[0] + nchrobs[1] + nchrobs[2] + nchrobs[3] + nchrobs[4];
                double nchrobsEUR = nchrobs[5] + nchrobs[6] + nchrobs[7] + nchrobs[8] + nchrobs[9];
                double nchrobsAFR = nchrobs[10] + nchrobs[11] + nchrobs[12] + nchrobs[13] + nchrobs[14] + nchrobs[15] + nchrobs[16];
                double nchrobsAMR = nchrobs[17] + nchrobs[18] + nchrobs[19] + nchrobs[20];
                double nchrobsSAS = nchrobs[21] + nchrobs[22] + nchrobs[23] + nchrobs[24] + nchrobs[25];
                
                double mafEAS = macEAS/nchrobsEAS;
                double mafEUR = macEUR/nchrobsEUR;
                double mafAFR = macAFR/nchrobsAFR;
                double mafAMR = macAMR/nchrobsAMR;
                double mafSAS = macSAS/nchrobsSAS;
                
                 double mafAll = (macEAS + macEUR + macAFR + macAMR + macSAS) / (nchrobsEAS + nchrobsEUR + nchrobsAFR + nchrobsAMR + nchrobsSAS);
                if(mafAll < mafCutoff){
                    break;
                }
                
                string outline = currentChr + "\t" + currentSNP + "\t" + currentA1 + "\t" + currentA2 + "\t";
                string outline2 = outline;
                for(size_t i=0; i<(maf.size()-1); i++){
                    for(size_t j=(i+1); j<maf.size(); j++){
                        outline = outline + to_string(maf[i] - maf[j]) + "\t";
                    }
                }
                
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
                numSNP++;
            }
            break;
        }
        
        vector<string> vsline = split(fline, " ", true);
        if(currentTag != vsline[1] + vsline[3] + vsline[4]){
            if(currentSNP != ""){
                if(numRace != population.size()){
                    cout << "Race information is missing for " + currentSNP <<endl;
                    currentChr = vsline[0];
                    currentSNP = vsline[1];
                    currentA1 = vsline[3];
                    currentA2 = vsline[4];
                    currentTag = currentSNP + currentA1 + currentA2;
                    numRace = 0;
                    continue;
                }
                
                double macEAS = mac[0] + mac[1] + mac[2] + mac[3] + mac[4];
                double macEUR = mac[5] + mac[6] + mac[7] + mac[8] + mac[9];
                double macAFR = mac[10] + mac[11] + mac[12] + mac[13] + mac[14] + mac[15] + mac[16];
                double macAMR = mac[17] + mac[18] + mac[19] + mac[20];
                double macSAS = mac[21] + mac[22] + mac[23] + mac[24] + mac[25];
                
                double nchrobsEAS = nchrobs[0] + nchrobs[1] + nchrobs[2] + nchrobs[3] + nchrobs[4];
                double nchrobsEUR = nchrobs[5] + nchrobs[6] + nchrobs[7] + nchrobs[8] + nchrobs[9];
                double nchrobsAFR = nchrobs[10] + nchrobs[11] + nchrobs[12] + nchrobs[13] + nchrobs[14] + nchrobs[15] + nchrobs[16];
                double nchrobsAMR = nchrobs[17] + nchrobs[18] + nchrobs[19] + nchrobs[20];
                double nchrobsSAS = nchrobs[21] + nchrobs[22] + nchrobs[23] + nchrobs[24] + nchrobs[25];
                
                double mafEAS = macEAS/nchrobsEAS;
                double mafEUR = macEUR/nchrobsEUR;
                double mafAFR = macAFR/nchrobsAFR;
                double mafAMR = macAMR/nchrobsAMR;
                double mafSAS = macSAS/nchrobsSAS;
                
                double mafAll = (macEAS + macEUR + macAFR + macAMR + macSAS) / (nchrobsEAS + nchrobsEUR + nchrobsAFR + nchrobsAMR + nchrobsSAS);
                if(mafAll < mafCutoff){
                    currentChr = vsline[0];
                    currentSNP = vsline[1];
                    currentA1 = vsline[3];
                    currentA2 = vsline[4];
                    currentTag = currentSNP + currentA1 + currentA2;
                    numRace = 0;
                    
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
                    continue;
                }
                
                mafMajor[0] = mafEAS;
                mafMajor[1] = mafEUR;
                mafMajor[2] = mafAFR;
                mafMajor[3] = mafAMR;
                mafMajor[4] = mafSAS;
                
                string outline = currentChr + "\t" + currentSNP + "\t" + currentA1 + "\t" + currentA2 + "\t";
                string outline2 = outline;
                for(size_t i=0; i<(maf.size()-1); i++){
                    for(size_t j=(i+1); j<maf.size(); j++){
                        outline = outline + to_string(maf[i] - maf[j]) + "\t";
                    }
                }
                
                
                
                for(size_t i=0; i<(mafMajor.size()-1); i++){
                    for(size_t j=(i+1); j<mafMajor.size(); j++){
                        outline2 = outline2 + to_string(mafMajor[i] - mafMajor[j]) + "\t";
                    }
                }
                
                fout1 << outline << endl;
                fout2 << outline2 << endl;
                numSNP++;
                if(numSNP % 10000 == 0){
                    cout << numSNP << endl;
                }
            }
            currentChr = vsline[0];
            currentSNP = vsline[1];
            currentA1 = vsline[3];
            currentA2 = vsline[4];
            currentTag = currentSNP + currentA1 + currentA2;
            numRace = 0;
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
    
    cout << "Total SNP: " << numSNP << endl;
    
    return true;
}


bool processAllele(const char* alleleFile, const char* snpInfoFile, const char* outputFile1, const char* outputFile2, double mafCutoff) {
    
    // SNP informaiton
    vector<int> rsSNPInfo;
    vector<string> SNPInfo;
    
    ifstream fin(snpInfoFile);
    if(!fin.is_open()){
        cout << "Cannot open file " << snpInfoFile << endl;
        return false;
    }
    
    string header;
    getline(fin, header);
    
    while(!fin.eof()){
        string fline;
        getline(fin, fline);
        if(fline.size() < 2){
            break;
        }
        
        size_t idx = fline.find("\t");
        int rsTmp = stoi(fline.substr(0, idx));
        string snpInfoTmp = fline.substr(idx+1);
        rsSNPInfo.push_back(rsTmp);
        SNPInfo.push_back(snpInfoTmp);
    }
    fin.close();
    fin.clear();
    
    
    
    // alele information
    fin.open(alleleFile);
    if(!fin.is_open()){
        cerr << "Cannot ope file "<<alleleFile<<endl;
        return false;
    }
    
    // inlcude position, gene and other information
    ofstream fout1(outputFile1);
    ofstream fout2(outputFile2);
    
    fout1 << "CHR\tPos\tSNP\tA1\tA2\tGene\tInfo\t";
    for(size_t i = 0; i<(population.size()-1); i++){
        for(size_t j=(i+1); j<population.size(); j++){
            fout1 << population[i] << "-" << population[j] << "\t";
        }
    }
    fout1 <<endl;
    
    fout2 << "CHR\tPos\tSNP\tA1\tA2\tGene\tInfo\t";
    for(size_t i = 0; i<(majorPopulation.size()-1); i++){
        for(size_t j=(i+1); j<majorPopulation.size(); j++){
            fout2 << majorPopulation[i] << "-" << majorPopulation[j] << "\t";
        }
    }
    fout2 <<endl;
    
    vector<double> maf(population.size(), 0.0);
    //maf of 5 major population
    vector<double> mafMajor(5, 0.0);
    size_t numRace = 0;
    vector<double> mac(population.size(), 0.0);
    vector<double> nchrobs(population.size(), 0.0);
    string currentChr = "";
    string currentPos = "";
    string currentSNP = "";
    string currentA1 = "";
    string currentA2 = "";
    string currentGene = "";
    string currentInfo = "";
    string currentTag = "";
    
    
    getline(fin, header); // remove header
    size_t numSNP = 0;
    size_t numMiss = 0;
    while (!fin.eof()) {
        string fline;
        getline(fin, fline);
        
        if(fline.size() < 2){
            // process the last record
            if(numRace != population.size()){
                cout << "Race information is missing for " + currentSNP <<endl;
            } else {
                double macEAS = mac[0] + mac[1] + mac[2] + mac[3] + mac[4];
                double macEUR = mac[5] + mac[6] + mac[7] + mac[8] + mac[9];
                double macAFR = mac[10] + mac[11] + mac[12] + mac[13] + mac[14] + mac[15] + mac[16];
                double macAMR = mac[17] + mac[18] + mac[19] + mac[20];
                double macSAS = mac[21] + mac[22] + mac[23] + mac[24] + mac[25];
                
                double nchrobsEAS = nchrobs[0] + nchrobs[1] + nchrobs[2] + nchrobs[3] + nchrobs[4];
                double nchrobsEUR = nchrobs[5] + nchrobs[6] + nchrobs[7] + nchrobs[8] + nchrobs[9];
                double nchrobsAFR = nchrobs[10] + nchrobs[11] + nchrobs[12] + nchrobs[13] + nchrobs[14] + nchrobs[15] + nchrobs[16];
                double nchrobsAMR = nchrobs[17] + nchrobs[18] + nchrobs[19] + nchrobs[20];
                double nchrobsSAS = nchrobs[21] + nchrobs[22] + nchrobs[23] + nchrobs[24] + nchrobs[25];
                
                double mafEAS = macEAS/nchrobsEAS;
                double mafEUR = macEUR/nchrobsEUR;
                double mafAFR = macAFR/nchrobsAFR;
                double mafAMR = macAMR/nchrobsAMR;
                double mafSAS = macSAS/nchrobsSAS;
                
                 double mafAll = (macEAS + macEUR + macAFR + macAMR + macSAS) / (nchrobsEAS + nchrobsEUR + nchrobsAFR + nchrobsAMR + nchrobsSAS);
                if(mafAll < mafCutoff){
                    break;
                }
                
                string outline = currentChr + "\t" + currentPos + "\t" + currentSNP + "\t" + currentA1 + "\t" + currentA2 + "\t" + currentGene + "\t" + currentInfo + "\t";
                string outline2 = outline;
                for(size_t i=0; i<(maf.size()-1); i++){
                    for(size_t j=(i+1); j<maf.size(); j++){
                        outline = outline + to_string(maf[i] - maf[j]) + "\t";
                    }
                }
                
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
                if(currentPos == ""){
                    numMiss++;
                }
                numSNP++;
            }
            break;
        }
        
        vector<string> vsline = split(fline, " ", true);
        if(currentTag != vsline[1] + vsline[3] + vsline[4]){
            if(currentSNP != ""){
                if(numRace != population.size()){
                    cout << "Race information is missing for " + currentSNP <<endl;
                    currentChr = vsline[0];
                    currentSNP = vsline[1];
                    currentA1 = vsline[3];
                    currentA2 = vsline[4];
                    currentTag = currentSNP + currentA1 + currentA2;
                    
                    int rsNumber = stoi(currentSNP.substr(2));
                    int idx = binSearch(rsSNPInfo, rsNumber);
                    
                    if(idx<0){
                        currentPos = "";
                        currentGene = "";
                        currentInfo = "";
                    } else {
                        vector<string> vsTmp = split(SNPInfo[idx], "\t");
                        currentPos = vsTmp[0];
                        currentGene = vsTmp[1];
                        currentInfo = vsTmp[2];
                    }
                    
                    numRace = 0;
                    continue;
                }
                
                double macEAS = mac[0] + mac[1] + mac[2] + mac[3] + mac[4];
                double macEUR = mac[5] + mac[6] + mac[7] + mac[8] + mac[9];
                double macAFR = mac[10] + mac[11] + mac[12] + mac[13] + mac[14] + mac[15] + mac[16];
                double macAMR = mac[17] + mac[18] + mac[19] + mac[20];
                double macSAS = mac[21] + mac[22] + mac[23] + mac[24] + mac[25];
                
                double nchrobsEAS = nchrobs[0] + nchrobs[1] + nchrobs[2] + nchrobs[3] + nchrobs[4];
                double nchrobsEUR = nchrobs[5] + nchrobs[6] + nchrobs[7] + nchrobs[8] + nchrobs[9];
                double nchrobsAFR = nchrobs[10] + nchrobs[11] + nchrobs[12] + nchrobs[13] + nchrobs[14] + nchrobs[15] + nchrobs[16];
                double nchrobsAMR = nchrobs[17] + nchrobs[18] + nchrobs[19] + nchrobs[20];
                double nchrobsSAS = nchrobs[21] + nchrobs[22] + nchrobs[23] + nchrobs[24] + nchrobs[25];
                
                double mafEAS = macEAS/nchrobsEAS;
                double mafEUR = macEUR/nchrobsEUR;
                double mafAFR = macAFR/nchrobsAFR;
                double mafAMR = macAMR/nchrobsAMR;
                double mafSAS = macSAS/nchrobsSAS;
                
                double mafAll = (macEAS + macEUR + macAFR + macAMR + macSAS) / (nchrobsEAS + nchrobsEUR + nchrobsAFR + nchrobsAMR + nchrobsSAS);
                if(mafAll < mafCutoff){
                    currentChr = vsline[0];
                    currentSNP = vsline[1];
                    currentA1 = vsline[3];
                    currentA2 = vsline[4];
                    currentTag = currentSNP + currentA1 + currentA2;
                    
                    int rsNumber = stoi(currentSNP.substr(2));
                    int idx = binSearch(rsSNPInfo, rsNumber);
                    
                    if(idx<0){
                        currentPos = "";
                        currentGene = "";
                        currentInfo = "";
                    } else {
                        vector<string> vsTmp = split(SNPInfo[idx], "\t");
                        currentPos = vsTmp[0];
                        currentGene = vsTmp[1];
                        currentInfo = vsTmp[2];
                    }
                    
                    numRace = 0;
                    
                    idx = -1;
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
                    continue;
                }
                
                mafMajor[0] = mafEAS;
                mafMajor[1] = mafEUR;
                mafMajor[2] = mafAFR;
                mafMajor[3] = mafAMR;
                mafMajor[4] = mafSAS;
                
                string outline = currentChr + "\t" + currentPos + "\t" + currentSNP + "\t" + currentA1 + "\t" + currentA2 + "\t" + currentGene + "\t" + currentInfo + "\t";
                string outline2 = outline;
                for(size_t i=0; i<(maf.size()-1); i++){
                    for(size_t j=(i+1); j<maf.size(); j++){
                        outline = outline + to_string(maf[i] - maf[j]) + "\t";
                    }
                }
                
                
                
                for(size_t i=0; i<(mafMajor.size()-1); i++){
                    for(size_t j=(i+1); j<mafMajor.size(); j++){
                        outline2 = outline2 + to_string(mafMajor[i] - mafMajor[j]) + "\t";
                    }
                }
                
                fout1 << outline << endl;
                fout2 << outline2 << endl;
                if(currentPos == ""){
                    numMiss++;
                }
                numSNP++;
                if(numSNP % 10000 == 0){
                    cout << numSNP << endl;
                }
            }
            currentChr = vsline[0];
            currentSNP = vsline[1];
            currentA1 = vsline[3];
            currentA2 = vsline[4];
            currentTag = currentSNP + currentA1 + currentA2;
            
            int rsNumber = stoi(currentSNP.substr(2));
            int idx = binSearch(rsSNPInfo, rsNumber);
            
            if(idx<0){
                currentPos = "";
                currentGene = "";
                currentInfo = "";
            } else {
                vector<string> vsTmp = split(SNPInfo[idx], "\t");
                currentPos = vsTmp[0];
                currentGene = vsTmp[1];
                currentInfo = vsTmp[2];
            }
            
            numRace = 0;
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
    
    cout << "SNP cannot find: " << numMiss << endl;
    cout << "Total SNP: " << numSNP << endl;
    
    return true;
}

