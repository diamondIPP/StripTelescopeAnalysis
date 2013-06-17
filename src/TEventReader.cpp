/*
 * TEventReader.cpp
 *
 *  Created on: 10.11.2011
 *      Author: bachmair
 */

#include "../include/TEventReader.hh"

TEventReader::TEventReader(std::string pathName) {
	verbosity=0;
	current_event = 0;
	store_threshold = 2;
	tree =NULL;
	file=NULL;

	SetTree(pathName);//tree);
	initialiseTree();
	if(!this->isOK()){
		cout<<"TADCEventReader::TADCEventReader is not correctly initialized.. EXIT PROGRAM"<<endl;
		exit(-1);
	}
	if (verbosity) cout<<"PedTree Entries():"<<PedTree->GetEntries()<<endl;
	this->GetEvent(0);
}

TEventReader::~TEventReader() {
	// TODO Auto-generated destructor stub
}


void TEventReader::setTree(){

}
