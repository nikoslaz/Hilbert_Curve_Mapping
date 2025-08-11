#ifndef MESH_H
#define MESH_H

#include <cstdint>
#include <iostream>
#include <vector>
#include <utility> 
#include <fstream> 
#include <tuple>   
#include <algorithm> 
#include <cmath>     
#include <stdexcept> 
#include <string>    

using namespace std;

class Mesh {

public:
    int dimensions;
    int bits;
    vector<vector<uint64_t>> mesh; 
    vector<int> num_each_dim;      

    Mesh(vector<int> dimensions_vec) {
        if (dimensions_vec.empty()) {
            throw std::runtime_error("Mesh constructor requires non-empty dimensions vector.");
        }
        this->num_each_dim = dimensions_vec;
        this->dimensions = dimensions_vec.size();
        // Calculate bits needed based on the largest dimension size
        int max_dim_size = *max_element(num_each_dim.begin(), num_each_dim.end());
        if (max_dim_size <= 0) {
             throw std::runtime_error("Mesh dimensions must be positive.");
        }
        // Calculate bits needed to represent indices up to max_dim_size - 1
        this->bits = (max_dim_size == 1) ? 1 : static_cast<int>(ceil(log2(static_cast<double>(max_dim_size))));
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------- //
    
    /* Recursively generates all points in a grid defined by pointsPerDimension */
    void generate_points_recursive(vector<uint64_t>& current, int dim, const vector<int>& pointsPerDimension, vector<vector<uint64_t>>& mesh_output) {
        if (dim == current.size()) {
            // Base case: Store the generated point
            mesh_output.push_back(current);
            return;
        }

        if (dim < 0 || dim >= pointsPerDimension.size()) {
             throw std::out_of_range("Dimension index out of range in generate_points_recursive.");
        }

        // Iterate over the number of points allowed in this dimension
        for (int i = 0; i < pointsPerDimension[dim]; ++i) {
            current[dim] = static_cast<uint64_t>(i);
            generate_points_recursive(current, dim + 1, pointsPerDimension, mesh_output);
        }
    }


    void generate_points() {
        if (dimensions <= 0 || num_each_dim.empty()) {
            cerr << "Warning: Cannot generate points with invalid dimensions." << endl;
            return;
        }
        mesh.clear(); // Clear previous points
        vector<uint64_t> current_point(dimensions); // Temporary vector for recursion
        generate_points_recursive(current_point, 0, num_each_dim, mesh);
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------- //

    void print_points() {
        cout << "\nGenerated " << mesh.size() << " points in " << dimensions << "D space:\n";
        for (const auto& point : mesh) {
            cout << "(";
            for (size_t i = 0; i < point.size(); ++i) {
                cout << point[i] << (i + 1 < point.size() ? ", " : "");
            }
            cout << ")\n";
        }
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------- //

    vector<pair<pair<int, int>, pair<int, int>>> generate_vertexes(int X, int Y) {
        if (X <= 0 || Y <= 0) {
             throw std::runtime_error("generate_vertexes(2D): Dimensions X and Y must be positive.");
        }
        vector<pair<pair<int, int>, pair<int, int>>> result;
        // Reserve space (approximate, slightly less than X*Y*2 without wrap-around)
        result.reserve((X - 1) * Y + X * (Y - 1));

        for (int i = 0; i < X; i++) {
            for (int j = 0; j < Y; j++) {
                // Horizontal edge (only if not at the right boundary)
                if (j + 1 < Y) {
                    result.push_back(make_pair(make_pair(i, j), make_pair(i, j + 1)));
                }

                // Vertical edge (only if not at the bottom boundary)
                if (i + 1 < X) {
                    result.push_back(make_pair(make_pair(i, j), make_pair(i + 1, j)));
                }
            }
        }
        return result;
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------- //

    vector<pair<tuple<int, int, int>, tuple<int, int, int>>> generate_vertexes(int X, int Y, int Z) {
         if (X <= 0 || Y <= 0 || Z <= 0) {
             throw std::runtime_error("generate_vertexes(3D): Dimensions X, Y, and Z must be positive.");
        }
        vector<pair<tuple<int, int, int>, tuple<int, int, int>>> result;
        // Reserve space (approximate, slightly less than X*Y*Z*3 without wrap-around)
        result.reserve((X - 1) * Y * Z + X * (Y - 1) * Z + X * Y * (Z - 1));

        for (int i = 0; i < X; i++) {
            for (int j = 0; j < Y; j++) {
                for (int k = 0; k < Z; k++) {
                    // X-direction edge (only if not at the X boundary)
                    if (i + 1 < X) {
                        result.push_back(make_pair(make_tuple(i, j, k), make_tuple(i + 1, j, k)));
                    }

                    // Y-direction edge (only if not at the Y boundary)
                    if (j + 1 < Y) {
                        result.push_back(make_pair(make_tuple(i, j, k), make_tuple(i, j + 1, k)));
                    }

                    // Z-direction edge (only if not at the Z boundary)
                    if (k + 1 < Z) {
                        result.push_back(make_pair(make_tuple(i, j, k), make_tuple(i, j, k + 1)));
                    }
                }
            }
        }
        return result;
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------- //

    void export_mesh_csv(const vector<pair<pair<int, int>, pair<int, int>>>& vertexes, const string& filename = "bin/vertexes_2d.csv") {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file for writing: " << filename << endl;
            return;
        }

        file << "x1,y1,x2,y2\n"; // Header

        for (const auto& edge : vertexes) {
            file << edge.first.first << "," << edge.first.second << ","
                 << edge.second.first << "," << edge.second.second << "\n";
        }

        file.close();
        cout << "2D mesh edges saved to " << filename << endl;
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------- //
    
    void export_mesh_csv(const vector<pair<tuple<int, int, int>, tuple<int, int, int>>>& vertexes, const string& filename = "bin/vertexes_3d.csv") {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file for writing: " << filename << endl;
            return;
        }

        file << "x1,y1,z1,x2,y2,z2\n"; // Header

        for (const auto& edge : vertexes) {
            const auto& first_tuple = edge.first;
            const auto& second_tuple = edge.second;

            file << get<0>(first_tuple) << "," << get<1>(first_tuple) << "," << get<2>(first_tuple) << ","
                 << get<0>(second_tuple) << "," << get<1>(second_tuple) << "," << get<2>(second_tuple) << "\n";
        }

        file.close();
        cout << "3D mesh edges saved to " << filename << endl;
    }
};

#endif //MESH_H