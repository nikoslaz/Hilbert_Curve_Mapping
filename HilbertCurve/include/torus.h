#include <iostream>
#include <vector>
#include <utility>
#include <tuple>
#include <fstream>

using namespace std;

class Torus {
private:
    int X, Y, Z; // Dimensions
    bool is3D;   // Flag to distinguish 2D vs 3D
    vector<vector<pair<int, int>>> adjList2D;       // 2D adjacency list
    vector<vector<tuple<int, int, int>>> adjList3D; // 3D adjacency list


    // Convert 2D coordinates to 1D index
    int toIndex2D(int i, int j) const {
        return i * Y + j;
    }

    // Convert 3D coordinates to 1D index
    int toIndex3D(int i, int j, int k) const {
        return i * Y * Z + j * Z + k;
    }

    // Normalize coordinates with wrap-around
    int mod(int x, int size) const {
        return (x % size + size) % size; // Handles negative numbers
    }

public:
    vector<vector<uint64_t>> points;

    // 2D Constructor
    Torus(int X_, int Y_) : X(X_), Y(Y_), Z(0), is3D(false) {
        points.resize(X * Y);
        for (int i = 0; i < X; i++) {
            for (int j = 0; j < Y; j++) {
                int current = toIndex2D(i, j);
                points[current] = {
                    static_cast<uint64_t>(toIndex2D(mod(i + 1, X), j)), // Right
                    static_cast<uint64_t>(toIndex2D(mod(i - 1, X), j)), // Left
                    static_cast<uint64_t>(toIndex2D(i, mod(j + 1, Y))), // Up
                    static_cast<uint64_t>(toIndex2D(i, mod(j - 1, Y)))  // Down
                };
            }
        }
    }

    // 3D Constructor
    Torus(int X_, int Y_, int Z_) : X(X_), Y(Y_), Z(Z_), is3D(true) {
        points.resize(X * Y * Z);
        for (int i = 0; i < X; i++) {
            for (int j = 0; j < Y; j++) {
                for (int k = 0; k < Z; k++) {
                    int current = toIndex3D(i, j, k);
                    points[current] = {
                        static_cast<uint64_t>(toIndex3D(mod(i + 1, X), j, k)), // Right
                        static_cast<uint64_t>(toIndex3D(mod(i - 1, X), j, k)), // Left
                        static_cast<uint64_t>(toIndex3D(i, mod(j + 1, Y), k)), // Up
                        static_cast<uint64_t>(toIndex3D(i, mod(j - 1, Y), k)), // Down
                        static_cast<uint64_t>(toIndex3D(i, j, mod(k + 1, Z))), // Forward
                        static_cast<uint64_t>(toIndex3D(i, j, mod(k - 1, Z)))  // Backward
                    };
                }
            }
        }
    }

    void generate_points(int dimensions, int X, int Y, int Z, int n_nodes) {
        vector<vector<uint64_t>> points(n_nodes);
        int idx = 0;
        if (dimensions == 2) {
            if (Z == 1) {
                // 2D grid in X-Y plane
                for (int i = 0; i < X && idx < n_nodes; i++) {
                    for (int j = 0; j < Y && idx < n_nodes; j++) {
                        points[idx] = {static_cast<uint64_t>(i), static_cast<uint64_t>(j)};
                        idx++;
                    }
                }
            } else if (Y == 1) {
                // 2D grid in X-Z plane
                for (int i = 0; i < X && idx < n_nodes; i++) {
                    for (int k = 0; k < Z && idx < n_nodes; k++) {
                        points[idx] = {static_cast<uint64_t>(i), static_cast<uint64_t>(k)};
                        idx++;
                    }
                }
            } else { // X == 1
                // 2D grid in Y-Z plane
                for (int j = 0; j < Y && idx < n_nodes; j++) {
                    for (int k = 0; k < Z && idx < n_nodes; k++) {
                        points[idx] = {static_cast<uint64_t>(j), static_cast<uint64_t>(k)};
                        idx++;
                    }
                }
            }
        } else {
            // 3D grid
            for (int i = 0; i < X && idx < n_nodes; i++) {
                for (int j = 0; j < Y && idx < n_nodes; j++) {
                    for (int k = 0; k < Z && idx < n_nodes; k++) {
                        points[idx] = {static_cast<uint64_t>(i), static_cast<uint64_t>(j), static_cast<uint64_t>(k)};
                        idx++;
                    }
                }
            }
        }
    }

    // Get neighbors for 2D torus
    vector<pair<int, int>> getNeighbors2D(int i, int j) const {
        if (!is3D && i >= 0 && i < X && j >= 0 && j < Y) {
            return adjList2D[toIndex2D(i, j)];
        }
        return {};
    }

    // Get neighbors for 3D torus
    vector<tuple<int, int, int>> getNeighbors3D(int i, int j, int k) const {
        if (is3D && i >= 0 && i < X && j >= 0 && j < Y && k >= 0 && k < Z) {
            return adjList3D[toIndex3D(i, j, k)];
        }
        return {};
    }

    // Print the torus structure
    void printTorus() const {
        if (!is3D) {
            // 2D Torus
            cout << "2D Torus (" << X << " x " << Y << "):\n";
            cout << "Index | Coordinates | Neighbors (Right, Left, Up, Down)\n";
            cout << "--------------------------------------------\n";
            for (int i = 0; i < X; i++) {
                for (int j = 0; j < Y; j++) {
                    int idx = toIndex2D(i, j);
                    const auto& neighbors = points[idx];
                    cout << idx << " | (" << i << ", " << j << ") | ("
                              << neighbors[0] << ", " << neighbors[1] << ", "
                              << neighbors[2] << ", " << neighbors[3] << ")\n";
                }
            }
        } else {
            // 3D Torus
            cout << "3D Torus (" << X << " x " << Y << " x " << Z << "):\n";
            cout << "Index | Coordinates | Neighbors (Right, Left, Up, Down, Forward, Backward)\n";
            cout << "------------------------------------------------------------\n";
            for (int i = 0; i < X; i++) {
                for (int j = 0; j < Y; j++) {
                    for (int k = 0; k < Z; k++) {
                        int idx = toIndex3D(i, j, k);
                        const auto& neighbors = points[idx];
                        cout << idx << " | (" << i << ", " << j << ", " << k << ") | ("
                                  << neighbors[0] << ", " << neighbors[1] << ", "
                                  << neighbors[2] << ", " << neighbors[3] << ", "
                                  << neighbors[4] << ", " << neighbors[5] << ")\n";
                    }
                }
            }
        }
    }

    void export2DToCSV(const string& filename) const {
        ofstream file(filename);
        file << "x1,y1,x2,y2\n";
        // Export Right and Up edges to avoid duplicates
        for (int i = 0; i < X; i++) {
            for (int j = 0; j < Y; j++) {
                // Right edge: (i, j) to (i+1, j)
                int rightI = mod(i + 1, X);
                file << i << "," << j << "," << rightI << "," << j << "\n";
                // Up edge: (i, j) to (i, j+1)
                int upJ = mod(j + 1, Y);
                file << i << "," << j << "," << i << "," << upJ << "\n";
            }
        }
        file.close();
    }
    
    void export3DToCSV(const string& filename) const {
        ofstream file(filename);
        file << "x1,y1,z1,x2,y2,z2\n";
        for (int i = 0; i < X; i++) {
            for (int j = 0; j < Y; j++) {
                for (int k = 0; k < Z; k++) {
                    // Right edge: (i, j, k) to (i+1, j, k)
                    int rightI = mod(i + 1, X);
                    file << i << "," << j << "," << k << "," << rightI << "," << j << "," << k << "\n";
                    // Up edge: (i, j, k) to (i, j+1, k)
                    int upJ = mod(j + 1, Y);
                    file << i << "," << j << "," << k << "," << i << "," << upJ << "," << k << "\n";
                    // Forward edge: (i, j, k) to (i, j, k+1)
                    int forwardK = mod(k + 1, Z);
                    file << i << "," << j << "," << k << "," << i << "," << j << "," << forwardK << "\n";
                }
            }
        }
        file.close();
    }
                    
};