#include <torus.h>
#include <torus_donut.h>
#include <hilbert_curve.h>
#include <dist.h>

int main(int argc, char* argv[]) {
    // Parse Arguments
    AppConfig config = parseArguments(argc, argv);
    cout << "Configuration: X=" << config.X << ", Y=" << config.Y << ", Z=" << config.Z
         << ", Dim=" << config.dimensions << ", File=" << config.trafficFilePath << endl;

    // Read Traffic Data
    int n_nodes = 0;
    vector<TrafficEntry> entries = readTrafficData(config.trafficFilePath, n_nodes);
    cout << "Read " << entries.size() << " traffic entries. Detected " << n_nodes << " nodes." << endl;

    // Validate Node Count
    // if (config.n_nodes != n_nodes) {
    //     throw runtime_error("Error: Dimensions X*Y*Z (" + to_string(config.n_nodes) +
    //     ") do not match the number of nodes found in CSV (" + to_string(n_nodes) + ")");
    // }

    // Initialize Original Torus and Export Points
    // Torus torus(config.dimensions == 2 ? (config.Z == 1 ? Torus(config.X, config.Y) : (config.Y == 1 ? Torus(config.X, config.Z) : Torus(config.Y, config.Z))) : Torus(config.X, config.Y, config.Z));
    // torus.export3DToCSV("bin/locs_torus.csv");
    // // Initialize points array
    // torus.generate_points(config.dimensions, config.X, config.Y, config.Z, config.n_nodes);

    // Generate Hilbert Curve Mapping
    HilbertData hilbertData = generateHilbertMapping(config);
    cout << "Generated Hilbert mapping for " << hilbertData.sortedPoints.size() << " points." << endl;

    // Export Hilbert Points
    Hilbert_Curve hilbert_exporter;
    hilbert_exporter.export_hilbert_csv(hilbertData.sortedPoints);

    // Analyze Traffic Locality
    // cout << "\nTorus traffic comparison (Before vs. After Hilbert Mapping):\n";
    // analyzeHighTrafficPairsTorus(entries, hilbertData.newCoords, torus.points, config.X, config.Y, config.Z, 0);
    
    // Torus Donut Visualization 
    float majorRadius = 5.0f;
    float minorRadius = 3.0f;

    Torus_Donut visualTorus(config.X, config.Y, config.Z, majorRadius, minorRadius);

    int surface_u_resolution = 60; 
    int surface_v_resolution = 30; 
    visualTorus.generateFullTorusSurface(surface_u_resolution, surface_v_resolution);
    visualTorus.saveFullTorusSurfaceToCSV("torus_surface.csv");

    visualTorus.clearTransformedPoints(); 
    int steps_per_hilbert_segment = 5;
    visualTorus.mapHilbertPoints(hilbertData, steps_per_hilbert_segment); 
    visualTorus.savePointsToCSV("torus_points.csv"); 

    return 0;
}
