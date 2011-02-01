/*
 * main.cpp
 * Diamond
 *
 * Created by Lukas Baeni on 19.01.11.
 */

//#include "TROOT.h"
//#include "Clustering.class.cpp"
#include "SlidingPedestal.class.cpp"
#include "Clustering.class.cpp"
#include <fstream>
#include <iostream>

using namespace std;

int main () {
	cout << "hello world!" << endl;
	SlidingPedestal sl(15208);
	sl.Slide(10000);
	
	Clustering cl(15208);
	cl.ClusterRun(1, 0);
}

//void readSettings() {
	
//}