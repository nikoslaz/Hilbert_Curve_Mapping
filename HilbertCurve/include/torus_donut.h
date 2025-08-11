#ifndef TORUS_DONUT_H
#define TORUS_DONUT_H

#include <structs.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm> 

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

typedef struct { 
    float x, y, z;
} Point3D;

class Torus_Donut {
private:
    int gridX = 0, gridY = 0, gridZ = 0;
    float map_outerRad = 0.0f;
    float map_innerRad = 0.0f;

    vector<Point3D> transformed_points;
    vector<Point3D> full_torus_surface_points;
    int surface_u_steps = 0;
    int surface_v_steps = 0;

public:
    Torus_Donut(int gx, int gy, int gz, float targetOuterRad, float targetInnerRad)
    : gridX(gx), gridY(gy), gridZ(gz), map_outerRad(targetOuterRad), map_innerRad(targetInnerRad) {
        if (this->gridX <= 0) this->gridX = 1;
        if (this->gridY <= 0) this->gridY = 1;
        if (this->gridZ <= 0) this->gridZ = 1; 
    }

    void clearTransformedPoints() {
        transformed_points.clear();
    }

    float parametric_curve(float val0, float val1, float t) {
        return val0 * (1.0f - t) + val1 * t;
    }

    void mapHilbertPoints(HilbertData& hilbertData, int steps_per_segment) {
        const float N = static_cast<float>(gridX);
        const float n = static_cast<float>(gridY);

        for (size_t i = 0; i < hilbertData.sortedPoints.size() - 1; ++i) {
            const auto& c_coords = hilbertData.sortedPoints[i].second;
            const auto& n_coords = hilbertData.sortedPoints[i+1].second;

            if (c_coords.empty() || n_coords.empty() || c_coords.size() < 2 || n_coords.size() < 2) {
                continue; 
            }

            float u_0 = (2.0f * M_PI *static_cast<float>(c_coords[0]))/N;
            float v_0 = (2.0f * M_PI *static_cast<float>(c_coords[1]))/n;

            float u_1 = (2.0f * M_PI *static_cast<float>(n_coords[0]))/N;
            float v_1 = (2.0f * M_PI *static_cast<float>(n_coords[1]))/n;

            if(u_1 - u_0 > M_PI) u_0+=2*M_PI;
            if(u_1 - u_0 < -M_PI) u_0-=2*M_PI;

            if(v_1 - v_0 > M_PI) v_0+=2*M_PI;
            if(v_1 - v_0 < -M_PI) v_0-=2*M_PI;

            for (int step = 0; step < steps_per_segment; ++step) {
                float t = static_cast<float>(step) / static_cast<float>(steps_per_segment);
                float u = parametric_curve(u_0, u_1, t);
                float v = parametric_curve(v_0, v_1, t);

                Point3D p;
                p.x = (map_outerRad + map_innerRad * cos(v)) * cos(u);
                p.y = (map_outerRad + map_innerRad * cos(v)) * sin(u); 
                p.z = map_innerRad * sin(v);
                transformed_points.push_back(p);
            }

        }
        
        if (!hilbertData.sortedPoints.empty()) {
             const auto& l_coords = hilbertData.sortedPoints.back().second;
             if (!l_coords.empty() && l_coords.size() >=2) {
                float i = static_cast<float>(l_coords[0]);
                float j = static_cast<float>(l_coords[1]);

                float u = (2.0f * M_PI * i) / N;
                float v = (2.0f * M_PI * j) / n;
                Point3D p;
                p.x = (map_outerRad + map_innerRad * cos(v)) * cos(u);
                p.y = (map_outerRad + map_innerRad * cos(v)) * sin(u);
                p.z = map_innerRad * sin(v);
                transformed_points.push_back(p);
             }
        }
    }

    void generateFullTorusSurface(int u_res, int v_res) {
        full_torus_surface_points.clear();
        surface_u_steps = std::max(2, u_res); 
        surface_v_steps = std::max(2, v_res); 

        for (int i = 0; i < surface_u_steps; ++i) {
            float u = static_cast<float>(i) / static_cast<float>(surface_u_steps) * 2.0f * M_PI;
            for (int j = 0; j < surface_v_steps; ++j) {
                float v = static_cast<float>(j) / static_cast<float>(surface_v_steps) * 2.0f * M_PI;
                Point3D p;
                p.x = (map_outerRad + map_innerRad * cos(v)) * cos(u); 
                p.y = (map_outerRad + map_innerRad * cos(v)) * sin(u); 
                p.z = map_innerRad * sin(v);
                full_torus_surface_points.push_back(p);
            }
        }
    }
    
    void savePointsToCSV(const string& filename) const {
        ofstream outFile(filename);
        if (!outFile.is_open()) { return; }
        outFile << fixed << setprecision(6);
        outFile << "# Mapped Hilbert Data on Torus CSV (Nodes or Interpolated Path)\n";
        outFile << "# gridX=" << gridX << "\n";
        outFile << "# gridY=" << gridY << "\n";
        outFile << "# gridZ=" << gridZ << "\n";
        outFile << "# map_outerRad=" << map_outerRad << "\n";
        outFile << "# map_innerRad=" << map_innerRad << "\n";
        outFile << "x,y,z\n";
        for (const auto& p : transformed_points) {
            outFile << p.x << "," << p.y << "," << p.z << "\n";
        }
        outFile.close();
    }

    void saveFullTorusSurfaceToCSV(const string& filename) const {
        ofstream outFile(filename);
        if (!outFile.is_open()) { return; }
        outFile << fixed << setprecision(6);
        outFile << "# Full Torus Surface Points CSV (for Visualization)\n";
        outFile << "# map_outerRad=" << map_outerRad << "\n";
        outFile << "# map_innerRad=" << map_innerRad << "\n";
        outFile << "# surface_u_steps=" << surface_u_steps << "\n";
        outFile << "# surface_v_steps=" << surface_v_steps << "\n";
        outFile << "x,y,z\n";
        for (const auto& p : full_torus_surface_points) {
            outFile << p.x << "," << p.y << "," << p.z << "\n";
        }
        outFile.close();
    }
};

#endif // TORUS_DONUT_H