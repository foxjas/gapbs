#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <functional>
#include <climits>
#include <ctime>
#include <stdlib.h>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <sys/time.h>

#include "GorderGraph.h"
#include "GorderUtil.h"

using namespace std;

const int INPUTNUM=1;

int main(int argc, char* argv[]){
	ios::sync_with_stdio(false);
	int i;
	int W=5;
    int k = 0;
    bool undirected=false;
	clock_t start, end;
	string filename;

	if(argc==1){
		cout << "please provide parameter" << endl;
		exit(0);
	}

	i=1;
	while(i<argc){
		if(strcmp("-w", argv[i])==0){
			i++;
			W=atoi(argv[i]);
			if(W<=0){
				cout << "w should be larger than 0" << endl;
				quit();
			}
			i++;
		}
        else if(strcmp("-k", argv[i])==0){
			i++;
			k=atoi(argv[i]);
			// if(k<=0){
			// 	cout << "k should be larger than 0" << endl;
			// 	quit();
			// }
			i++;
		}
        else if(strcmp("-u", argv[i])==0){
		    undirected=true;
			i++;
		}
		else{
			filename=argv[i++];
		}
	}

	srand(time(0));

	GorderGraph g, g_tmp;
	string name;
	name=extractFilename(filename.c_str());
	g.setFilename(name);

	//start=clock();
	g.readGorderGraph(filename, false); // want to output relabeled graph as undirected
	vector<int> removed;
	if (k) {
		g.RemoveGreaterThanTopK(0); // this is just to have consistent graph creation
		g_tmp.setFilename(name);
		g_tmp.readGorderGraph(filename, undirected);
    	removed = g_tmp.RemoveGreaterThanTopK(k);
	} else {
	    g.Transform();
	}

	cout << name << " readGorderGraph is complete." << endl;
	//end=clock();
	//cout << "Time Cost: " << (double)(end-start)/CLOCKS_PER_SEC << endl;

	start=clock();
	vector<int> order;
	if (k) {
		g_tmp.GorderGreedy(order, removed, W);
	} else {
		g.GorderGreedy(order, W);
	}
	end=clock();
	// cout << "New ordering size: " << order.size() << endl;
	cout << "ReOrdered Time Cost: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
	cout << "Begin Output the Reordered GorderGraph" << endl;
	g.PrintReOrderedGorderGraph(order);
	cout << endl;

}

