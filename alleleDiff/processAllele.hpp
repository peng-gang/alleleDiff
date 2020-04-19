//
//  processAllele.hpp
//  alleleDiff
//
//  Created by Gang Peng on 4/16/20.
//  Copyright © 2020 Gang Peng. All rights reserved.
//

#ifndef processAllele_hpp
#define processAllele_hpp

#include <stdio.h>

bool processAllele(const char* alleleFile, const char* outputFile1, const char* outputFile2, double mafCutoff = 0.01);

bool processAllele(const char* alleleFile, const char* snpInfoFile, const char* outputFile1, const char* outputFile2, double mafCutoff = 0.01);

#endif /* processAllele_hpp */
