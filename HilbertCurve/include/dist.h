#ifndef DISTANCE_H
#define DISTANCE_H

#include <structs.h>

#define MANH_AFTER(dx, dy, dz) ((dx) + (dy) + (dz))
#define EUCL_AFTER(dx, dy, dz) (sqrt((dx) * (dx) + (dy) * (dy) + (dz) * (dz)))
#define CHEB_AFTER(dx, dy, dz) (max({(dx), (dy), (dz)}))

using namespace std;

// ------------------------------------------------------------------------------------------------------------------------------------------- //

AppConfig parseArguments(int argc, char* argv[]) {
    if (argc < 5) {
        throw runtime_error("Usage: " + string(argv[0]) + " <X> <Y> <Z> <traffic_csv_filename>");
    }

    AppConfig config;
    config.X = atoi(argv[1]);
    config.Y = atoi(argv[2]);
    config.Z = atoi(argv[3]);
    config.trafficFilePath = string("../../test_files/") + argv[4]; 

    if (config.X <= 0 || config.Y <= 0 || config.Z <= 0) {
        throw runtime_error("Dimensions X, Y, Z must be positive integers.");
    }

    config.dimensions = (config.X == 1 || config.Y == 1 || config.Z == 1) ? 2 : 3;
    config.n_nodes = config.X * config.Y * config.Z;

    cout << "Input Args: X=" << argv[1] << ", Y=" << argv[2] << ", Z=" << argv[3] << ", File=" << argv[4] << endl;

    return config;
}

// ------------------------------------------------------------------------------------------------------------------------------------------- //

vector<TrafficEntry> readTrafficData(const string& filePath, int& n_nodes) {
    ifstream file(filePath);
    if (!file.is_open()) {
        throw runtime_error("Cannot open traffic CSV file: " + filePath);
    }

    vector<TrafficEntry> entries;
    string line;
    getline(file, line); // Skip header

    int max_node = -1;
    int lineNum = 1;
    while (getline(file, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        int src, dst;
        double traffic;
        char comma;
        if (!(ss >> src >> comma >> dst >> comma >> traffic)) {
            cerr << "Error parsing line: " << line << endl;
            continue;
        }
        entries.push_back({src, dst, traffic});
        max_node = max(max_node, max(src, dst));
    }

    if (max_node == -1) {
        n_nodes = 0; 
    } else {
        n_nodes = max_node + 1;
    }

    return entries;
}

// ------------------------------------------------------------------------------------------------------------------------------------------- //

bool compare_hilbert_points(const pair<uint64_t, vector<uint64_t>>& a, const pair<uint64_t, vector<uint64_t>>& b) {
    return a.first < b.first;
}

// ------------------------------------------------------------------------------------------------------------------------------------------- //

HilbertData generateHilbertMapping(const AppConfig& config) {
    HilbertData data;
    Hilbert_Curve hilbert; // Instance for encoding

    int max_dim_val = max({config.X, config.Y, config.Z});
    int bits = 0;
    while ((1 << bits) < max_dim_val) bits++;

    data.newCoords.resize(config.n_nodes);
    vector<vector<uint64_t>> original_coords(config.n_nodes);

    int index = 0;
    vector<vector<uint64_t>> current_point_coords(1, vector<uint64_t>(config.dimensions)); // Pre-allocate for encode

    if (config.dimensions == 2) {
        // 2D handling - assumes Z=1 means XY plane, Y=1 means XZ, X=1 means YZ
        int dim1 = (config.Z == 1 || config.Y == 1) ? config.X : config.Y;
        int dim2 = (config.Z == 1) ? config.Y : config.Z; // If Z=1 use Y, else use Z

        for (int i = 0; i < dim1 && index < config.n_nodes; ++i) {
            for (int j = 0; j < dim2 && index < config.n_nodes; ++j) {
                vector<uint64_t> coords_2d = {static_cast<uint64_t>(i), static_cast<uint64_t>(j)};
                current_point_coords[0] = coords_2d;
                uint64_t h = hilbert.encode(current_point_coords, config.dimensions, bits)[0];

                vector<uint64_t> coords_3d;
                if (config.Z == 1) coords_3d = {coords_2d[0], coords_2d[1], 0};         // XY Plane
                else if (config.Y == 1) coords_3d = {coords_2d[0], 0, coords_2d[1]};    // XZ Plane
                else coords_3d = {0, coords_2d[0], coords_2d[1]};                       // YZ Plane

                data.sortedPoints.push_back({h, coords_3d});
                index++;
            }
        }
    } else { // 3D
        for (int i = 0; i < config.X && index < config.n_nodes; ++i) {
            for (int j = 0; j < config.Y && index < config.n_nodes; ++j) {
                for (int k = 0; k < config.Z && index < config.n_nodes; ++k) {
                    vector<uint64_t> coords_3d = {static_cast<uint64_t>(i), static_cast<uint64_t>(j), static_cast<uint64_t>(k)};
                    current_point_coords[0] = coords_3d;
                    uint64_t h = hilbert.encode(current_point_coords, config.dimensions, bits)[0];

                    data.sortedPoints.push_back({h, coords_3d});
                    index++;
                }
            }
        }
    }

    // Sort by Hilbert index
    sort(data.sortedPoints.begin(), data.sortedPoints.end(), compare_hilbert_points);

    // Create the mapping from Hilbert index (position in sorted list) to the coordinate
    for (int h_idx = 0; h_idx < config.n_nodes; ++h_idx) {
        if (h_idx < data.sortedPoints.size()) {
             data.newCoords[h_idx] = data.sortedPoints[h_idx].second;
        } else {
            // Handle potential mismatch if index didn't reach n_nodes
             cerr << "Warning: Hilbert points generated (" << data.sortedPoints.size()
                  << ") mismatch expected nodes (" << config.n_nodes << "). Index reached: " << index << endl;
             // Assign a default or handle error appropriately
             if (!data.newCoords[h_idx].empty()) data.newCoords[h_idx].assign(config.dimensions, 0);
        }
    }

    return data;
}

// ------------------------------------------------------------------------------------------------------------------------------------------- //

double euclideanDistance(const vector<uint64_t>& a, const vector<uint64_t>& b) {
    if (a.size() != b.size()) {
        cerr << "Error: Vectors must have the same size for Euclidean distance." << endl;
        return -1.0;
    }
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        double diff = static_cast<double>(a[i]) - static_cast<double>(b[i]);
        sum += diff * diff;
    }
    return sqrt(sum);
}

double euclideanDistanceTorus(int src_idx, int dst_idx, int X, int Y, int Z) {
    int src_x = src_idx / (Y * Z);
    int src_y = (src_idx / Z) % Y;
    int src_z = src_idx % Z;
    int dst_x = dst_idx / (Y * Z);
    int dst_y = (dst_idx / Z) % Y;
    int dst_z = dst_idx % Z;

    int dx = min(abs(src_x - dst_x), X - abs(src_x - dst_x));
    int dy = min(abs(src_y - dst_y), Y - abs(src_y - dst_y));
    int dz = min(abs(src_z - dst_z), Z - abs(src_z - dst_z));

    return sqrt(static_cast<double>(dx * dx + dy * dy + dz * dz));
}

// ------------------------------------------------------------------------------------------------------------------------------------------- //

double manhattanDistance(const vector<uint64_t>& a, const vector<uint64_t>& b) {
    if (a.size() != b.size()) {
        cerr << "Error: Vectors must have the same size for Manhattan distance." << endl;
        return -1.0;
    }
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        sum += abs(static_cast<double>(a[i]) - static_cast<double>(b[i]));
    }
    return sum;
}

double manhattanDistanceTorus(int src_idx, int dst_idx, int X, int Y, int Z) {
    int src_x = src_idx / (Y * Z);
    int src_y = (src_idx / Z) % Y;
    int src_z = src_idx % Z;
    int dst_x = dst_idx / (Y * Z);
    int dst_y = (dst_idx / Z) % Y;
    int dst_z = dst_idx % Z;

    int dx = min(abs(src_x - dst_x), X - abs(src_x - dst_x));
    int dy = min(abs(src_y - dst_y), Y - abs(src_y - dst_y));
    int dz = min(abs(src_z - dst_z), Z - abs(src_z - dst_z));

    return static_cast<double>(dx + dy + dz);
}

// ------------------------------------------------------------------------------------------------------------------------------------------- //

double chebyshevDistance(const vector<uint64_t>& a, const vector<uint64_t>& b) {
    if (a.size() != b.size()) {
        cerr << "Error: Vectors must have the same size for Chebyshev distance." << endl;
        return -1.0;
    }
    double maxDiff = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        double diff = abs(static_cast<double>(a[i]) - static_cast<double>(b[i]));
        maxDiff = max(maxDiff, diff);
    }
    return maxDiff;
}

double chebyshevDistanceTorus(int src_idx, int dst_idx, int X, int Y, int Z) {
    int src_x = src_idx / (Y * Z);
    int src_y = (src_idx / Z) % Y;
    int src_z = src_idx % Z;
    int dst_x = dst_idx / (Y * Z);
    int dst_y = (dst_idx / Z) % Y;
    int dst_z = dst_idx % Z;

    int dx = min(abs(src_x - dst_x), X - abs(src_x - dst_x));
    int dy = min(abs(src_y - dst_y), Y - abs(src_y - dst_y));
    int dz = min(abs(src_z - dst_z), Z - abs(src_z - dst_z));

    return static_cast<double>(max({dx, dy, dz}));
}

// ------------------------------------------------------------------------------------------------------------------------------------------- //

void export_dist_csv(vector<double> eucl, vector<double> manh, vector<double> cheb) {
    ofstream file("bin/dists.csv");
    file << "eucl, manh, cheb" << endl;
    for (int i = 0; i < eucl.size(); i++) {
        file << eucl[i] << ", " << manh[i] << ", " << cheb[i] << endl;
    }
}

// ------------------------------------------------------------------------------------------------------------------------------------------- //

void compareTrafficBeforeAndAfter(const vector<TrafficEntry>& entries, 
    const vector<int>& hilbert_order_mapping, 
    const vector<pair<uint64_t, vector<uint64_t>>>& hilbertPoints, 
    int X, int Y, int Z) {
    // Before Hilbert Curve
    double total_traffic_before = 0.0;
    vector<double> eucl_before, manh_before, cheb_before;

    for (const auto& entry : entries) {
        total_traffic_before += entry.traffic;

        // 3D coordinates before Hilbert mapping
        int src_x = entry.src % X;
        int src_y = (entry.src / X) % Y;
        int src_z = entry.src / (X * Y);
        int dst_x = entry.dst % X;
        int dst_y = (entry.dst / X) % Y;
        int dst_z = entry.dst / (X * Y);

        vector<uint64_t> src_coords = {static_cast<uint64_t>(src_x), static_cast<uint64_t>(src_y), static_cast<uint64_t>(src_z)};
        vector<uint64_t> dst_coords = {static_cast<uint64_t>(dst_x), static_cast<uint64_t>(dst_y), static_cast<uint64_t>(dst_z)};
        eucl_before.push_back(euclideanDistance(src_coords, dst_coords));
        manh_before.push_back(manhattanDistance(src_coords, dst_coords));
        cheb_before.push_back(chebyshevDistance(src_coords, dst_coords));
    }

    double avg_eucl_before = accumulate(eucl_before.begin(), eucl_before.end(), 0.0) / eucl_before.size();
    double avg_manh_before = accumulate(manh_before.begin(), manh_before.end(), 0.0) / manh_before.size();
    double avg_cheb_before = accumulate(cheb_before.begin(), cheb_before.end(), 0.0) / cheb_before.size();

    // After Hilbert Curve
    double total_traffic_after = 0.0;
    vector<double> eucl_after, manh_after, cheb_after;

    for (const auto& entry : entries) {
        total_traffic_after += entry.traffic;

        int src_hilbert = hilbert_order_mapping[entry.src];
        int dst_hilbert = hilbert_order_mapping[entry.dst];

        auto src_it = find_if(hilbertPoints.begin(), hilbertPoints.end(),
            [src_hilbert](const auto& p) { return p.first == src_hilbert; });
        auto dst_it = find_if(hilbertPoints.begin(), hilbertPoints.end(),
            [dst_hilbert](const auto& p) { return p.first == dst_hilbert; });

        eucl_after.push_back(euclideanDistance(src_it->second, dst_it->second));
        manh_after.push_back(manhattanDistance(src_it->second, dst_it->second));
        cheb_after.push_back(chebyshevDistance(src_it->second, dst_it->second));
    }

    double avg_eucl_after = accumulate(eucl_after.begin(), eucl_after.end(), 0.0) / eucl_after.size();
    double avg_manh_after = accumulate(manh_after.begin(), manh_after.end(), 0.0) / manh_after.size();
    double avg_cheb_after = accumulate(cheb_after.begin(), cheb_after.end(), 0.0) / cheb_after.size();

    // Output results (aligned)
    cout << setw(45) << left << "Traffic Comparison:" << endl << endl;
    
    cout << setw(45) << left << "Before Hilbert Curve:" << endl;
    cout << setw(45) << left << "  Average Euclidean Distance:" << avg_eucl_before << endl;
    cout << setw(45) << left << "  Average Manhattan Distance:" << avg_manh_before << endl;
    cout << setw(45) << left << "  Average Chebyshev Distance:" << avg_cheb_before << endl << endl;

    cout << setw(45) << left << "After Hilbert Curve:" << endl;
    cout << setw(45) << left << "  Average Euclidean Distance:" << avg_eucl_after << endl;
    cout << setw(45) << left << "  Average Manhattan Distance:" << avg_manh_after << endl;
    cout << setw(45) << left << "  Average Chebyshev Distance:" << avg_cheb_after << endl << endl;

    int improved_eucl = 0;
    for (size_t i = 0; i < eucl_before.size(); ++i) {
        if (eucl_after[i] < eucl_before[i]) improved_eucl++;
    }
    cout << setw(55) << left 
         << "Number of traffic entries with improved locality (Euclidean):" 
         << improved_eucl << " out of " << eucl_before.size() << endl;

    int improved_manh = 0;
    for (size_t i = 0; i < manh_before.size(); ++i) {
        if (manh_after[i] < manh_before[i]) improved_manh++;
    }
    cout << setw(55) << left 
         << "Number of traffic entries with improved locality (Manhattan):" 
         << improved_manh << " out of " << manh_before.size() << endl;

    int improved_cheb = 0;
    for (size_t i = 0; i < cheb_before.size(); ++i) {
        if (cheb_after[i] < cheb_before[i]) improved_cheb++;
    }
    cout << setw(55) << left 
         << "Number of traffic entries with improved locality (Chebyshev):" 
         << improved_cheb << " out of " << cheb_before.size() << endl;
}

void compareTrafficBeforeAndAfterTorus(const vector<TrafficEntry>& entries, 
    const vector<int>& hilbert_order_mapping, 
    const vector<pair<uint64_t, vector<uint64_t>>>& hilbertPoints,
    const vector<vector<uint64_t>>& torusPoints,
    int X, int Y, int Z) {
    
    // "Before" using torusPoints (from the Torus object)
    double total_traffic_before = 0.0;
    vector<double> eucl_before, manh_before, cheb_before;

    for (const auto& entry : entries) {
        total_traffic_before += entry.traffic;

        eucl_before.push_back(euclideanDistanceTorus(entry.src, entry.dst, X, Y, Z));
        manh_before.push_back(manhattanDistanceTorus(entry.src, entry.dst, X, Y, Z));
        cheb_before.push_back(chebyshevDistanceTorus(entry.src, entry.dst, X, Y, Z));
    }

    double avg_eucl_before = accumulate(eucl_before.begin(), eucl_before.end(), 0.0) / eucl_before.size();
    double avg_manh_before = accumulate(manh_before.begin(), manh_before.end(), 0.0) / manh_before.size();
    double avg_cheb_before = accumulate(cheb_before.begin(), cheb_before.end(), 0.0) / cheb_before.size();

    // After Hilbert Curve, compute distances using the Hilbert-transformed points
    double total_traffic_after = 0.0;
    vector<double> eucl_after, manh_after, cheb_after;

    for (const auto& entry : entries) {
        total_traffic_after += entry.traffic;
    
        int src_hilbert = hilbert_order_mapping[entry.src];
        int dst_hilbert = hilbert_order_mapping[entry.dst];
    
        auto src_it = find_if(hilbertPoints.begin(), hilbertPoints.end(),
            [src_hilbert](const auto& p) { return p.first == src_hilbert; });
        auto dst_it = find_if(hilbertPoints.begin(), hilbertPoints.end(),
            [dst_hilbert](const auto& p) { return p.first == dst_hilbert; });
    
        eucl_after.push_back(euclideanDistance(src_it->second, dst_it->second));
        manh_after.push_back(manhattanDistance(src_it->second, dst_it->second));
        cheb_after.push_back(chebyshevDistance(src_it->second, dst_it->second));
    }

    double avg_eucl_after = accumulate(eucl_after.begin(), eucl_after.end(), 0.0) / eucl_after.size();
    double avg_manh_after = accumulate(manh_after.begin(), manh_after.end(), 0.0) / manh_after.size();
    double avg_cheb_after = accumulate(cheb_after.begin(), cheb_after.end(), 0.0) / cheb_after.size();

    // Aligning output using setw and left
    cout << setw(45) << left << "Traffic Comparison:" << endl << endl;
    
    cout << setw(45) << left << "Before Hilbert Curve (from Torus):" << endl;
    cout << setw(45) << left << "  Average Euclidean Distance:" << avg_eucl_before << endl;
    cout << setw(45) << left << "  Average Manhattan Distance:" << avg_manh_before << endl;
    cout << setw(45) << left << "  Average Chebyshev Distance:" << avg_cheb_before << endl << endl;
    
    cout << setw(45) << left << "After Hilbert Curve:" << endl;
    cout << setw(45) << left << "  Average Euclidean Distance:" << avg_eucl_after << endl;
    cout << setw(45) << left << "  Average Manhattan Distance:" << avg_manh_after << endl;
    cout << setw(45) << left << "  Average Chebyshev Distance:" << avg_cheb_after << endl << endl;
    
    int improved_eucl = 0;
    for (size_t i = 0; i < eucl_before.size(); ++i) {
        if (eucl_after[i] < eucl_before[i])
            improved_eucl++;
    }
    cout << setw(55) << left 
         << "Number of traffic entries with improved locality (Euclidean):" 
         << improved_eucl << " out of " << eucl_before.size() << endl;
    
    int improved_manh = 0;
    for (size_t i = 0; i < manh_before.size(); ++i) {
        if (manh_after[i] < manh_before[i])
            improved_manh++;
    }
    cout << setw(55) << left 
         << "Number of traffic entries with improved locality (Manhattan):" 
         << improved_manh << " out of " << manh_before.size() << endl;
    
    int improved_cheb = 0;
    for (size_t i = 0; i < cheb_before.size(); ++i) {
        if (cheb_after[i] < cheb_before[i])
            improved_cheb++;
    }
    cout << setw(55) << left 
         << "Number of traffic entries with improved locality (Chebyshev):" 
         << improved_cheb << " out of " << cheb_before.size() << endl;
}

// ------------------------------------------------------------------------------------------------------------------------------------------- //

void analyzeHighTrafficPairs(const vector<TrafficEntry>& entries, 
const vector<int>& hilbert_order_mapping, 
const vector<pair<uint64_t, vector<uint64_t>>>& hilbertPoints, 
int X, int Y, int Z, double traffic_threshold = 0.0) {
    map<pair<int, int>, double> pair_traffic;
    for (const auto& entry : entries) {
        pair<int, int> node_pair = {entry.src, entry.dst};
        pair_traffic[node_pair] += entry.traffic;
    }

    double total_eucl_byte_before = 0.0, total_eucl_byte_after = 0.0;
    double total_manh_byte_before = 0.0, total_manh_byte_after = 0.0;
    double total_cheb_byte_before = 0.0, total_cheb_byte_after = 0.0;
    double max_eucl_byte_before   = 0.0, max_eucl_byte_after = 0.0;
    double max_manh_byte_before   = 0.0, max_manh_byte_after = 0.0;
    double max_cheb_byte_before   = 0.0, max_cheb_byte_after = 0.0;

    cout << "High-Traffic Pair Analysis (Threshold: " << traffic_threshold << "):" << endl;

    for (const auto& [node_pair, traffic] : pair_traffic) {
        if (traffic < traffic_threshold) continue;

        int src = node_pair.first;
        int dst = node_pair.second;

        // Before Hilbert Curve (3D coordinates)
        int src_x = src % X;
        int src_y = (src / X) % Y;
        int src_z = src / (X * Y);
        int dst_x = dst % X;
        int dst_y = (dst / X) % Y;
        int dst_z = dst / (X * Y);

        vector<uint64_t> src_coords = {static_cast<uint64_t>(src_x), static_cast<uint64_t>(src_y), static_cast<uint64_t>(src_z)};
        vector<uint64_t> dst_coords = {static_cast<uint64_t>(dst_x), static_cast<uint64_t>(dst_y), static_cast<uint64_t>(dst_z)};

        double eucl_before = euclideanDistance(src_coords, dst_coords);
        double manh_before = manhattanDistance(src_coords, dst_coords);
        double cheb_before = chebyshevDistance(src_coords, dst_coords);

        double eucl_byte_before = eucl_before * traffic;
        double manh_byte_before = manh_before * traffic;
        double cheb_byte_before = cheb_before * traffic;

        total_eucl_byte_before += eucl_byte_before;
        total_manh_byte_before += manh_byte_before;
        total_cheb_byte_before += cheb_byte_before;
        max_eucl_byte_before = max(max_eucl_byte_before, eucl_byte_before);
        max_manh_byte_before = max(max_manh_byte_before, manh_byte_before);
        max_cheb_byte_before = max(max_cheb_byte_before, cheb_byte_before);

        // After Hilbert Curve
        int src_hilbert = hilbert_order_mapping[src];
        int dst_hilbert = hilbert_order_mapping[dst];
        auto src_it = find_if(hilbertPoints.begin(), hilbertPoints.end(), [src_hilbert](const auto& p) { return p.first == src_hilbert; });
        auto dst_it = find_if(hilbertPoints.begin(), hilbertPoints.end(), [dst_hilbert](const auto& p) { return p.first == dst_hilbert; });

        double eucl_after = -1.0, manh_after = -1.0, cheb_after = -1.0;
        double eucl_byte_after = -1.0, manh_byte_after = -1.0, cheb_byte_after = -1.0;

        if (src_it != hilbertPoints.end() && dst_it != hilbertPoints.end()) {
            eucl_after = euclideanDistance(src_it->second, dst_it->second);
            manh_after = manhattanDistance(src_it->second, dst_it->second);
            cheb_after = chebyshevDistance(src_it->second, dst_it->second);

            eucl_byte_after = eucl_after * traffic;
            manh_byte_after = manh_after * traffic;
            cheb_byte_after = cheb_after * traffic;

            total_eucl_byte_after += eucl_byte_after;
            total_manh_byte_after += manh_byte_after;
            total_cheb_byte_after += cheb_byte_after;
            max_eucl_byte_after = max(max_eucl_byte_after, eucl_byte_after);
            max_manh_byte_after = max(max_manh_byte_after, manh_byte_after);
            max_cheb_byte_after = max(max_cheb_byte_after, cheb_byte_after);
        } else {
            cerr << "Error: Node not found for src=" << src << " or dst=" << dst << endl;
        }
    }

    auto calc_percent = [](double before, double after) {
        return (before > 0) ? ((after - before) / before) * 100 : 0.0;
    };

    double total_eucl_percent = calc_percent(total_eucl_byte_before, total_eucl_byte_after);
    double total_manh_percent = calc_percent(total_manh_byte_before, total_manh_byte_after);
    double total_cheb_percent = calc_percent(total_cheb_byte_before, total_cheb_byte_after);
    double max_eucl_percent   = calc_percent(max_eucl_byte_before, max_eucl_byte_after);
    double max_manh_percent   = calc_percent(max_manh_byte_before, max_manh_byte_after);
    double max_cheb_percent   = calc_percent(max_cheb_byte_before, max_cheb_byte_after);

    cout << "\nDistance x Byte Summary:" << endl;
    cout << setw(35) << left << "Euclidean:" << endl;
    cout << setw(35) << left << "  Total Distance x Byte Before:" << fixed << setprecision(2) << total_eucl_byte_before << endl;
    cout << setw(35) << left << "  Total Distance x Byte After:"  << total_eucl_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Total):"    << total_eucl_percent << "% (" << (total_eucl_percent < 0 ? "Better" : "Worse") << ")" << endl;
    cout << setw(35) << left << "  Max Distance x Byte Before:"   << max_eucl_byte_before << endl;
    cout << setw(35) << left << "  Max Distance x Byte After:"    << max_eucl_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Max):"      << max_eucl_percent << "% (" << (max_eucl_percent < 0 ? "Better" : "Worse") << ")" << endl;

    cout << setw(35) << left << "\nManhattan:" << endl;
    cout << setw(35) << left << "  Total Distance x Byte Before:" << total_manh_byte_before << endl;
    cout << setw(35) << left << "  Total Distance x Byte After:"  << total_manh_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Total):"    << total_manh_percent << "% (" << (total_manh_percent < 0 ? "Better" : "Worse") << ")" << endl;
    cout << setw(35) << left << "  Max Distance x Byte Before:"   << max_manh_byte_before << endl;
    cout << setw(35) << left << "  Max Distance x Byte After:"    << max_manh_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Max):"      << max_manh_percent << "% (" << (max_manh_percent < 0 ? "Better" : "Worse") << ")" << endl;

    cout << setw(35) << left << "\nChebyshev:" << endl;
    cout << setw(35) << left << "  Total Distance x Byte Before:" << total_cheb_byte_before << endl;
    cout << setw(35) << left << "  Total Distance x Byte After:"  << total_cheb_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Total):"    << total_cheb_percent << "% (" << (total_cheb_percent < 0 ? "Better" : "Worse") << ")" << endl;
    cout << setw(35) << left << "  Max Distance x Byte Before:"   << max_cheb_byte_before << endl;
    cout << setw(35) << left << "  Max Distance x Byte After:"    << max_cheb_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Max):"      << max_cheb_percent << "% (" << (max_cheb_percent < 0 ? "Better" : "Worse") << ")" << endl;
}

void analyzeHighTrafficPairsTorus(const vector<TrafficEntry>& entries, 
    const vector<vector<uint64_t>>& new_coords, 
    const vector<vector<uint64_t>>& torusPoints,
    int X, int Y, int Z, double traffic_threshold = 0.0) {
    
    // First, aggregate traffic for each node pair.
    map<pair<int, int>, double> pair_traffic;
    for (const auto& entry : entries) {
        pair<int, int> node_pair = {entry.src, entry.dst};
        pair_traffic[node_pair] += entry.traffic;
    }

    double total_eucl_byte_before = 0.0, total_eucl_byte_after = 0.0;
    double total_manh_byte_before = 0.0, total_manh_byte_after = 0.0;
    double total_cheb_byte_before = 0.0, total_cheb_byte_after = 0.0;
    double max_eucl_byte_before   = 0.0, max_eucl_byte_after   = 0.0;
    double max_manh_byte_before   = 0.0, max_manh_byte_after   = 0.0;
    double max_cheb_byte_before   = 0.0, max_cheb_byte_after   = 0.0;

    // Debug print for each pair with aligned formatting
    for (const auto& [node_pair, traffic] : pair_traffic) {
        if (traffic < traffic_threshold) continue;

        int src = node_pair.first;
        int dst = node_pair.second;

        // Before Hilbert Curve: use torusPoints (pre-transformation)
        vector<uint64_t> src_coords = torusPoints[src];
        vector<uint64_t> dst_coords = torusPoints[dst];

        double eucl_before = euclideanDistanceTorus(src, dst, X, Y, Z);
        double manh_before = manhattanDistanceTorus(src, dst, X, Y, Z);
        double cheb_before = chebyshevDistanceTorus(src, dst, X, Y, Z);

        double eucl_byte_before = eucl_before * traffic;
        double manh_byte_before = manh_before * traffic;
        double cheb_byte_before = cheb_before * traffic;

        total_eucl_byte_before += eucl_byte_before;
        total_manh_byte_before += manh_byte_before;
        total_cheb_byte_before += cheb_byte_before;

        max_eucl_byte_before = max(max_eucl_byte_before, eucl_byte_before);
        max_manh_byte_before = max(max_manh_byte_before, manh_byte_before);
        max_cheb_byte_before = max(max_cheb_byte_before, cheb_byte_before);

        // After Hilbert Curve: use the transformed coordinates from new_coords
        const auto& src_new_coords = new_coords[src];
        const auto& dst_new_coords = new_coords[dst];

        // Compute distances using torus-aware functions
        double eucl_after = 0.0, manh_after = 0.0, cheb_after = 0.0;
        int dx = min(abs(static_cast<int>(src_new_coords[0]) - static_cast<int>(dst_new_coords[0])), 
                     X - abs(static_cast<int>(src_new_coords[0]) - static_cast<int>(dst_new_coords[0])));
        int dy = min(abs(static_cast<int>(src_new_coords[1]) - static_cast<int>(dst_new_coords[1])), 
                     Y - abs(static_cast<int>(src_new_coords[1]) - static_cast<int>(dst_new_coords[1])));
        int dz = min(abs(static_cast<int>(src_new_coords[2]) - static_cast<int>(dst_new_coords[2])), 
                     Z - abs(static_cast<int>(src_new_coords[2]) - static_cast<int>(dst_new_coords[2])));

        manh_after = MANH_AFTER(dx, dy, dz);
        eucl_after = EUCL_AFTER(dx, dy, dz);
        cheb_after = CHEB_AFTER(dx, dy, dz);

        double eucl_byte_after = eucl_after * traffic;
        double manh_byte_after = manh_after * traffic;
        double cheb_byte_after = cheb_after * traffic;

        total_eucl_byte_after += eucl_byte_after;
        total_manh_byte_after += manh_byte_after;
        total_cheb_byte_after += cheb_byte_after;
        max_eucl_byte_after = max(max_eucl_byte_after, eucl_byte_after);
        max_manh_byte_after = max(max_manh_byte_after, manh_byte_after);
        max_cheb_byte_after = max(max_cheb_byte_after, cheb_byte_after);
    }

    // Output the summary with aligned formatting
    cout << endl << setw(50) << left << "Distance x Byte Summary:" << endl;
    cout << string(50, '=') << endl;

    // Euclidean Section
    cout << setw(20) << left << "Euclidean:" << endl;
    cout << setw(35) << left << "  Total Distance x Byte Before:" 
         << fixed << setprecision(2) << total_eucl_byte_before << endl;
    cout << setw(35) << left << "  Total Distance x Byte After:" 
         << fixed << setprecision(2) << total_eucl_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Total):" 
         << fixed << setprecision(2) << (total_eucl_byte_after - total_eucl_byte_before) / total_eucl_byte_before * 100 
         << "%" << endl;
    cout << setw(35) << left << "  Max Distance x Byte Before:" 
         << fixed << setprecision(2) << max_eucl_byte_before << endl;
    cout << setw(35) << left << "  Max Distance x Byte After:" 
         << fixed << setprecision(2) << max_eucl_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Max):" 
         << fixed << setprecision(2) << (max_eucl_byte_after - max_eucl_byte_before) / max_eucl_byte_before * 100 
         << "%" << endl;

    // Manhattan Section
    cout << endl << setw(20) << left << "Manhattan:" << endl;
    cout << setw(35) << left << "  Total Distance x Byte Before:" 
         << fixed << setprecision(2) << total_manh_byte_before << endl;
    cout << setw(35) << left << "  Total Distance x Byte After:" 
         << fixed << setprecision(2) << total_manh_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Total):" 
         << fixed << setprecision(2) << (total_manh_byte_after - total_manh_byte_before) / total_manh_byte_before * 100 
         << "%" << endl;
    cout << setw(35) << left << "  Max Distance x Byte Before:" 
         << fixed << setprecision(2) << max_manh_byte_before << endl;
    cout << setw(35) << left << "  Max Distance x Byte After:" 
         << fixed << setprecision(2) << max_manh_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Max):" 
         << fixed << setprecision(2) << (max_manh_byte_after - max_manh_byte_before) / max_manh_byte_before * 100 
         << "%" << endl;

    // Chebyshev Section
    cout << endl << setw(20) << left << "Chebyshev:" << endl;
    cout << setw(35) << left << "  Total Distance x Byte Before:" 
         << fixed << setprecision(2) << total_cheb_byte_before << endl;
    cout << setw(35) << left << "  Total Distance x Byte After:" 
         << fixed << setprecision(2) << total_cheb_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Total):" 
         << fixed << setprecision(2) << (total_cheb_byte_after - total_cheb_byte_before) / total_cheb_byte_before * 100 
         << "%" << endl;
    cout << setw(35) << left << "  Max Distance x Byte Before:" 
         << fixed << setprecision(2) << max_cheb_byte_before << endl;
    cout << setw(35) << left << "  Max Distance x Byte After:" 
         << fixed << setprecision(2) << max_cheb_byte_after << endl;
    cout << setw(35) << left << "  Percentage Change (Max):" 
         << fixed << setprecision(2) << (max_cheb_byte_after - max_cheb_byte_before) / max_cheb_byte_before * 100 
         << "%" << endl;
}

#endif