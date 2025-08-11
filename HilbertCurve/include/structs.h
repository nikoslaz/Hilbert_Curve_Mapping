#ifndef STRUCTS_H
#define STRUCTS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <map>
#include <iomanip>

using namespace std;

// Structure to hold traffic
struct TrafficEntry {
    int src;
    int dst;
    double traffic;   
};

// Structure to hold application configuration/arguments
struct AppConfig {
    int X, Y, Z;
    int dimensions;
    int n_nodes;
    string trafficFilePath;
};

// Structure to hold results from Hilbert generation
struct HilbertData {
    vector<pair<uint64_t, vector<uint64_t>>> sortedPoints;
    vector<vector<uint64_t>> newCoords;
};

#endif // STRUCTS_H