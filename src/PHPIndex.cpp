/*
 * PHPIndex.cpp
 *
 *  Created on: Sep 24, 2013
 *      Author: bachmair
 */

#include "../include/PHPIndex.hh"

PHPIndex::PHPIndex() {
    // TODO Auto-generated constructor stub
    this->subdirPath="";
    this->mainPath="";
    this->path="";
    setTitle("index");
    setFileName("index.php");
}

PHPIndex::~PHPIndex() {
    // TODO Auto-generated destructor stub
}



void PHPIndex::generateHTMLFile(){
    stringstream htmlOutputFileName;
//    if(verbosity>2)cout<<"FILE: "<<fileName<<endl;
    htmlOutputFileName<<fileGenPath<<"/"<<  fileName;
//    if(verbosity)cout<<"create HTML file: \""<<htmlOutputFileName.str()<<"\""<<endl;
    html_summary.open(htmlOutputFileName.str().c_str());
    fillIndex.
    if (verbosity>2)cout<<"::GenerateHTML()::close html_summary"<<flush;
    html_summary.close();
    if(verbosity>2)cout<<"...DONE"<<endl;
}
