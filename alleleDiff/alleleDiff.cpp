//
//  alleleDiff.cpp
//  alleleDiff
//
//  Created by Gang Peng on 4/15/20.
//  Copyright © 2020 Gang Peng. All rights reserved.
//

#include <iostream>
#include <fstream>

#include "normal.h"
#include "readFile.hpp"
#include "processAllele.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    if(argc < 5){
        cout << "Calculate allele frequency difference between race ethnicity groups" << endl;
        cout << "Usage 1: " << endl;
        cout << "alleleDiff  AlleleFile  OutputFile1 OutputFile2 mafCutoff" << endl;
        cout << "AlleleFile includes MAF for 26 populations" << endl;
        cout << "OutputFile1 is for allele frequency difference between 26 populations" << endl;
        cout << "OutputFile2 is for allele frequency difference betwwen 5 major populations" << endl;
        cout << "Usage 2: " << endl;
        cout << "alleleDiff  AlleleFile SNPInfoFile OutputFile1 OutputFile2 mafCutoff" << endl;
        cout << "AlleleFile includes MAF for 26 populations" << endl;
        cout << "SNPInfoFile includes SNP informaiton, such as location and associated genes" << endl;
        cout << "OutputFile1 is for allele frequency difference between 26 populations" << endl;
        cout << "OutputFile2 is for allele frequency difference betwwen 5 major populations" << endl;
        cout << "At least 4 papameters and less than 6 parameters should be set." << endl;
        return 0;
    } else if(argc == 5) {
        cout << "Start to calculate without SNP information" << endl;
        if(!processAllele(argv[1], argv[2], argv[3], stof(string(argv[4])))) {
            cout<<"Error during process" << endl;
            return 1;
        } else {
            cout<<"Finished" <<endl;
        }
    } else if(argc == 6) {
        cout << "Start to calculate with SNP information" << endl;
    } else {
        cout << "Too many parameters!" << endl;
        cout << "Usage 1: " << endl;
        cout << "alleleDiff  AlleleFile  OutputFile1 OutputFile2 mafCutoff" << endl;
        cout << "AlleleFile includes MAF for 26 populations" << endl;
        cout << "OutputFile1 is for allele frequency difference between 26 populations" << endl;
        cout << "OutputFile2 is for allele frequency difference betwwen 5 major populations" << endl;
        cout << "Usage 2: " << endl;
        cout << "alleleDiff  AlleleFile SNPInfoFile OutputFile1 OutputFile2 mafCutoff" << endl;
        cout << "AlleleFile includes MAF for 26 populations" << endl;
        cout << "SNPInfoFile includes SNP informaiton, such as location and associated genes" << endl;
        cout << "OutputFile1 is for allele frequency difference between 26 populations" << endl;
        cout << "OutputFile2 is for allele frequency difference betwwen 5 major populations" << endl;
        cout << "At least 4 papameters and less than 6 parameters should be set." << endl;
        return 0;
    }
    
    
    
    return 0;
}
