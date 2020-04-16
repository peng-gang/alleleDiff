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

bool processAllele(const char* alleleFile, const char* outputFile){
    ifstream fin(alleleFile);
    if(!fin.is_open()){
        cerr << "Cannot ope file "<<alleleFile<<endl;
        return false;
    }
    
    while (!fin.eof()) {
        string fline;
        string SNP = "";
        getline(fin, fline); // remove header
        getline(fin, fline);
        
        if(fline.size() < 2){
            break;
        }
        
        vector<string> vsline = split(fline, "\t");
        
    }
    return true;
}
