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
#include "populationInfo.h"
#include "readFile.hpp"
#include "processAllele.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    if(argc < 4){
        cout << "Calculate allele frequency difference between race ethnicity groups" << endl;
        cout << "Usage: " << endl;
        cout << "alleleDiff  AlleleFile SNPInfoFile OutputFile" << endl;
        cout << "All three parameters should be set." << endl;
        return 0;
    }
    
    
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}
