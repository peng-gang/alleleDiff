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
        cout << "Usage: " << endl;
        cout << "alleleDiff  AlleleFile SNPInfoFile OutputFile1 OutputFile2" << endl;
        cout << "OutputFile1 is for allele frequency difference between 26 populations" << endl;
        cout << "OutputFile2 is for allele frequency difference betwwen 5 major populations" << endl;
        cout << "All four parameters should be set." << endl;
        return 0;
    }
    
    if(!processAllele(argv[1], argv[3], argv[4])) {
        cout<<"Error during process" << endl;
        return 1;
    } else {
        cout<<"Finished" <<endl;
    }
    
    return 0;
}
