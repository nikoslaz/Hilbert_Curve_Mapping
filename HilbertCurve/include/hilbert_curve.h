#ifndef HILBERT_CURVE_H
#define HILBERT_CURVE_H

#include <cstddef>     
#include <iostream>     
#include <vector>       
#include <utility>      
#include <fstream>      
#include <cstdint>      
#include <stdexcept>    
#include <algorithm>    
#include <cmath>        

using namespace std;

class Hilbert_Curve {

public:
    Hilbert_Curve() = default;

    // ------------------------------------------------------------------------------------------------------------------------------------------- //
  
    void export_hilbert_csv(const vector<pair<uint64_t, vector<uint64_t>>>& hilbertPoints) const {
        ofstream file("locs_hilbert.csv");
        if (!file.is_open()) {
            cerr << "Error: Could not open file locs_hilbert.csv for writing." << endl;
            return;
        }

        if (hilbertPoints.empty()) {
            cerr << "Warning: No Hilbert points provided to export." << endl;
            file << "x,y,z\n"; // Write header even if empty
            file.close();
            return;
        }

        // Determine dimensions from the first point
        size_t num_dims = hilbertPoints[0].second.size();

        // Write header based on dimensions
        if (num_dims == 2) {
            file << "x,y\n";
        } else if (num_dims == 3) {
            file << "x,y,z\n";
        } else {
             // Write generic header for other dimensions
             for(size_t i = 0; i < num_dims; ++i) {
                 file << "dim" << i << (i == num_dims - 1 ? "" : ",");
             }
             file << "\n";
        }


        // Write points data
        for (const auto& p : hilbertPoints) {
            const auto& coords = p.second;
            if (coords.size() != num_dims) {
                 cerr << "Warning: Inconsistent dimensions found in hilbertPoints. Skipping point." << endl;
                 continue; // Skip inconsistent points
            }
            for (size_t i = 0; i < coords.size(); ++i) {
                file << coords[i] << (i == coords.size() - 1 ? "" : ",");
            }
            file << "\n";
        }
        file.close();
        cout << "Hilbert points saved to bin/locs_hilbert.csv" << endl;
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------- //
    // --- Gray Code and Bit Manipulation Helpers
    // ------------------------------------------------------------------------------------------------------------------------------------------- //

    /**
     * Converts a Gray-coded bit vector to a binary bit vector.
     * Rule: binary[0] = gray[0]; for i>=1, binary[i] = binary[i-1] XOR gray[i].
     */
    static vector<bool> gray2binary(const vector<bool>& gray) {
        vector<bool> binary(gray.size(), false);
        if (!gray.empty()) {
            binary[0] = gray[0];
            for (size_t i = 1; i < gray.size(); i++) {
                binary[i] = binary[i - 1] ^ gray[i];
            }
        }
        return binary;
    }

    /**
     * Applies gray2binary to each vector in a matrix.
     */
    static vector<vector<bool>> gray2binary(const vector<vector<bool>>& gray_matrix) {
        vector<vector<bool>> binary_matrix;
        binary_matrix.reserve(gray_matrix.size());
        for (const auto& gray_vec : gray_matrix) {
            binary_matrix.push_back(gray2binary(gray_vec));
        }
        return binary_matrix;
    }

    /**
     * Flattens a 3D Gray code matrix [N][num_dims][num_bits] into a 2D matrix [N][num_dims * num_bits].
     * This is often done after axis transposition in Hilbert encoding.
     */
    static vector<vector<bool>> flatten_gray(const vector<vector<vector<bool>>>& gray, int num_dims, int num_bits) {
        size_t N = gray.size();  // Number of locations
        if (N == 0 || gray[0].size() != num_dims || (num_dims > 0 && gray[0][0].size() != num_bits)) {
             throw runtime_error("Invalid dimensions for flatten_gray input.");
        }
        // Allocate a new 2D vector with shape [N][num_dims * num_bits]
        vector<vector<bool>> flatGray(N, vector<bool>(num_dims * num_bits, false));

        for (size_t i = 0; i < N; i++) {
            for (int b = 0; b < num_bits; b++) {
                for (int d = 0; d < num_dims; d++) {
                    // After swapping axes, the bit from the original gray[i][d][b]
                    // becomes the element at index (b * num_dims + d) in the flattened vector.
                    flatGray[i][b * num_dims + d] = gray[i][d][b];
                }
            }
        }
        return flatGray;
    }

    /**
     * Unpacks bytes into bits and truncates to the least significant num_bits.
     */
    static vector<bool> unpack_and_truncate(const vector<uint8_t>& bytes, int num_bits) {
        if (num_bits <= 0 || num_bits > bytes.size() * 8) {
            throw runtime_error("Invalid num_bits for unpack_and_truncate.");
        }
        // Unpack all bytes into bits.
        vector<bool> bits;
        bits.reserve(bytes.size() * 8);

        // Extract from most significant bit to least significant bit within each byte
        for (uint8_t byte : bytes) {
            for (int i = 7; i >= 0; --i) {
                bits.push_back((byte >> i) & 1);
            }
        }

        // Take only the last num_bits (least significant bits overall)
        vector<bool> truncated(bits.end() - num_bits, bits.end());
        return truncated;
    }

    /**
     * Packs a 64-bit vector (bool) into a uint64_t integer.
     * Assumes padded[0] is the most significant bit.
     */
    static uint64_t packBitsToUint64(const vector<bool>& padded) {
        if (padded.size() != 64) {
            throw runtime_error("packBitsToUint64 requires a vector with exactly 64 bits.");
        }

        uint64_t value = 0;
        // padded[0] is MSB, padded[63] is LSB
        for (int i = 0; i < 64; ++i) {
            if (padded[i]) {
                value |= (1ULL << (63 - i));
            }
        }
        return value;
    }

    /**
     * Applies packBitsToUint64 to each vector in a matrix.
     */
    static vector<uint64_t> convertAll(const vector<vector<bool>>& padded_matrix) {
        vector<uint64_t> result;
        result.reserve(padded_matrix.size());
        for (const auto& padded : padded_matrix) {
            result.push_back(packBitsToUint64(padded));
        }
        return result;
    }

    vector<uint64_t> encode(const vector<vector<uint64_t>>& locs, int num_dims, int num_bits) {
        // Keep the original shape for later
        int orig_shape = locs.size();
    
        for (const auto& loc : locs) {
            if (loc.size() != num_dims) {
                cerr << "Error. The shape of loc is not the same with number of dimensions. Must be the same!" << endl;
                exit(1);
            }
        }
    
        if (num_dims * num_bits > 64) {
            cerr << "Error.Total number of bits exceeds 64" << endl;
        }
    
        // Treat the location integers as 64-bit unsigned and then split them up into
        // a sequence of uint8s.  Preserve the association by dimension.
        // Allocate a 3D vector: (number of locations) x (num_dims) x (8 bytes per uint64_t)
        vector<vector<vector<uint8_t>>> locs_uint8(orig_shape, vector<vector<uint8_t>>(num_dims, vector<uint8_t>(8)));
        for (int i = 0; i < orig_shape; i++) {
            for (int j = 0; j < num_dims; j++) {
                for (int byte = 0; byte < 8; byte++) {
                    // Shift right so that the correct byte is in the lowest 8 bits.
                    locs_uint8[i][j][byte] = (uint8_t)((locs[i][j] >> (8 * (7 - byte))) & 0xFF);
                }
            }
        }
    
        // Now turn these into bits and truncate to num_bits.
        vector<vector<vector<bool>>> gray(orig_shape, vector<vector<bool>>(num_dims));
        for (int i = 0; i < orig_shape; i++) {
            for (int j = 0; j < num_dims; j++) {
                gray[i][j] = unpack_and_truncate(locs_uint8[i][j], num_bits);
            }
        }
    
        // Run the decoding process the other way.
        for(int bit = 0; bit < num_bits; bit++) {
            // Iterate forwards through the bits.
            for (int dim = 0; dim < num_dims; dim++) {
                // Iterate forwards through the dimensions.
                for (int idx = 0; idx < orig_shape; idx++) {
                    // Identify which ones have this bit active.
                    bool mask = gray[idx][dim][bit];
    
                    // Where this bit is on, invert the 0 dimension for lower bits.
                    if (mask) {
                        for (int b = bit + 1; b < num_bits; b++) {
                            // Flipping is equivalent to XOR with true.
                            gray[idx][0][b] = !gray[idx][0][b];
                        }
                    }
                    else {
                        for (int b = bit + 1; b < num_bits; b++) {
                            // Compute the XOR (difference) between dimension 0 and current dimension.
                            bool diff = gray[idx][0][b] ^ gray[idx][dim][b];
                            // If they differ, then this bit is marked to be flipped.
                            if (diff) {
                                // Flip both bits.
                                gray[idx][0][b]   = !gray[idx][0][b];
                                gray[idx][dim][b] = !gray[idx][dim][b];
                            }
                        }
                    }
                }
            }
        }
    
        // Now flatten out.
        vector<vector<bool>> flat_gray = flatten_gray(gray, num_dims, num_bits);
        
        //Convert Gray back to binary.
        vector<vector<bool>> hh_bin = gray2binary(flat_gray);
    
        // Pad back out to 64 bits.
        unsigned int extra_dims = 64 - num_bits * num_dims;
        vector<vector<bool>> padded;
        padded.reserve(hh_bin.size());
    
        for (const auto& row : hh_bin) {
            // Create a new padded row of length 64, initialized to false.
            vector<bool> paddedRow(64, false);
            // Copy the row into the padded row starting at index extra_dims.
            for (size_t i = 0; i < row.size(); i++) {
                paddedRow[extra_dims + i] = row[i];
            }
            padded.push_back(paddedRow);
        }
        
        vector<uint64_t> result = convertAll(padded);
    
        return result;
    }

};

#endif // HILBERT_CURVE_H