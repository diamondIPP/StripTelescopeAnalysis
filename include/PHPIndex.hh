/*
 * PHPIndex.hh
 *
 *  Created on: Sep 24, 2013
 *      Author: bachmair
 */

#ifndef PHPINDEX_HH_
#define PHPINDEX_HH_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TSystem.h"
#include "TDatime.h"
class PHPIndex {
public:
    PHPIndex();
    virtual ~PHPIndex();
    void setFileName(string Name);
    void setPathName(string pathName);
    void setMainPath(std::string mainPathName);
    void setSubdirPath(std::string subDirPath);
    void updatePath();
protected:
    string combineToString(string a,UInt_t b){
      stringstream output;
      output<<a<<b;
      return output.str();
    }
    string combineToString(string a,Double_t b,UInt_t precision=3){
      stringstream output;
      output<<a<<std::setprecision(precision)<<b;
      return;
    }
private:
    std::ofstream html_summary;
    TDatime dateandtime;
    std::string fileGenPath;
    std::string fileName;
    std::string path;
    std::string mainPath;
    std::string subdirPath;
};

#endif /* PHPINDEX_HH_ */
