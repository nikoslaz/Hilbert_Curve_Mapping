/*************************
 * TORUS DONUT CSV LOADING
 *************************/
const pathCsvFilename = '../../build/bin/torus_points.csv';
const surfaceCsvFilename = '../../build/bin/torus_surface.csv';
const plotDiv = document.getElementById('torus');


// --- Style Settings ---

const LINE_WIDTH = 8;
const SURFACE_OPACITY = 0.5;
const SURFACE_COLOR = 'rgb(200, 200, 200)';
const MARKER_SIZE = 1;
const OPACITY = 1;
const COLORSCALE = 'Viridis';

// --------------------------------

// Helper function to parse CSV text
function parseTorusCsv(csvText, filenameForError) {
    const data = {
        x: [], y: [], z: [],
        gridX: undefined, gridY: undefined, gridZ: undefined,
        surface_u_steps: undefined, surface_v_steps: undefined,
        map_outerRad: undefined, map_innerRad: undefined
    };
    const lines = csvText.trim().split('\n');
    let headerLineIndex = -1;

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line.startsWith('#')) {
            if (line.includes('gridX=')) data.gridX = parseInt(line.split('=')[1]);
            if (line.includes('gridY=')) data.gridY = parseInt(line.split('=')[1]);
            if (line.includes('gridZ=')) data.gridZ = parseInt(line.split('=')[1]);
            if (line.includes('surface_u_steps=')) data.surface_u_steps = parseInt(line.split('=')[1]);
            if (line.includes('surface_v_steps=')) data.surface_v_steps = parseInt(line.split('=')[1]);
            if (line.includes('map_outerRad=')) data.map_outerRad = parseFloat(line.split('=')[1]);
            if (line.includes('map_innerRad=')) data.map_innerRad = parseFloat(line.split('=')[1]);
        } else if (line.toLowerCase().startsWith('x,y,z')) {
            headerLineIndex = i;
            break;
        }
    }

    if (headerLineIndex === -1) {
        const errorMsg = `Could not find the 'x,y,z' header row in ${filenameForError}.`;
        console.error(errorMsg);
        throw new Error(errorMsg);
    }

    for (let i = headerLineIndex + 1; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line) {
            const values = line.split(',');
            if (values.length === 3) {
                data.x.push(parseFloat(values[0]));
                data.y.push(parseFloat(values[1]));
                data.z.push(parseFloat(values[2]));
            } else {
                console.warn(`  Skipping invalid data line ${i + 1} in ${filenameForError}: "${line}" (expected 3 values, got ${values.length})`);
            }
        }
    }
    
    return data;
}

Promise.all([
    fetch(surfaceCsvFilename).then(response => {
        if (!response.ok) throw new Error(`HTTP error! Status: ${response.status} - Could not fetch ${surfaceCsvFilename}`);
        return response.text();
    }).then(text => { console.log(`Successfully fetched ${surfaceCsvFilename}. Length: ${text.length}`); return text; }),
    fetch(pathCsvFilename).then(response => {
        if (!response.ok) throw new Error(`HTTP error! Status: ${response.status} - Could not fetch ${pathCsvFilename}`);
        return response.text();
    }).then(text => { console.log(`Successfully fetched ${pathCsvFilename}. Length: ${text.length}`); return text; })
])
.then(([surfaceCsvText, pathCsvText]) => { 
    const surfaceData = parseTorusCsv(surfaceCsvText, surfaceCsvFilename);
    const pathData = parseTorusCsv(pathCsvText, pathCsvFilename); 

    
    const allTraces = [];

    // --- Trace 1: Torus Surface ---
    if (surfaceData.x.length > 0 && surfaceData.surface_u_steps && surfaceData.surface_v_steps) {
        const u_steps = surfaceData.surface_u_steps;
        const v_steps = surfaceData.surface_v_steps;

        const x_surf = [];
        const y_surf = [];
        const z_surf = [];

        for (let i = 0; i < u_steps; i++) {
            x_surf.push(surfaceData.x.slice(i * v_steps, (i + 1) * v_steps));
            y_surf.push(surfaceData.y.slice(i * v_steps, (i + 1) * v_steps));
            z_surf.push(surfaceData.z.slice(i * v_steps, (i + 1) * v_steps));
        }
        
        const surfaceTrace = {
            x: x_surf, y: y_surf, z: z_surf,
            type: 'surface',
            colorscale: [[0, SURFACE_COLOR], [1, SURFACE_COLOR]],
            showscale: false,
            opacity: SURFACE_OPACITY,
            name: 'Torus Surface',
            showlegend: false,
            lighting: { ambient: 0.6, diffuse: 0.8, specular: 0.3, roughness: 0.4, fresnel: 0.2 },
            lightposition: { x: 100, y: 200, z: 50 }
        };
        allTraces.push(surfaceTrace);
    } else {
        console.warn("Could not generate Torus Surface trace. Missing data or u/v steps metadata.");
        if (!surfaceData.surface_u_steps) console.warn("  surface_u_steps is undefined.");
        if (!surfaceData.surface_v_steps) console.warn("  surface_v_steps is undefined.");
    }    

    // --- Trace 2: Mapped Hilbert Path Line ---
    if (pathData.x.length >= 2) {
        const mappedPathTrace = {
            x: pathData.x, y: pathData.y, z: pathData.z,
            mode: 'lines', type: 'scatter3d',
            line: {
                color: Array.from({length: pathData.x.length}, (_, k) => k),
                colorscale: COLORSCALE, showscale: true, width: LINE_WIDTH
            },
            name: 'Hilbert Path (on Torus)', showlegend: true
        };
        allTraces.push(mappedPathTrace);
    } else {
        console.warn("Not enough points to draw a Hilbert path line on the torus.");
    }
    
    if (allTraces.length === 0) {
        const errorMsg = "Error: No data to display for the torus plot. No traces were generated.";
        console.error(errorMsg);
        plotDiv.innerHTML = `<p>${errorMsg}</p>`;
        return;
    }

    const layout = {
        title: `Hilbert Path on Torus`,
        scene: {
            xaxis: { title: 'X Axis', showgrid: true, zeroline: true, backgroundcolor: "rgb(230, 230,230)" },
            yaxis: { title: 'Y Axis', showgrid: true, zeroline: true, backgroundcolor: "rgb(230, 230,230)" },
            zaxis: { title: 'Z Axis', showgrid: true, zeroline: true, backgroundcolor: "rgb(230, 230,230)" },
            aspectratio: { x: 1, y: 1, z: 0.8 },
            camera: { eye: {x: 1.8, y: 1.8, z: 0.8} }
        },
        margin: { l: 0, r: 0, b: 0, t: 40 },
        legend: { y: 0.9, x: 0.8 }
    };

    Plotly.newPlot(plotDiv, allTraces, layout);
})
.catch(error => {
    console.error('Error in Promise.all or subsequent .then for Torus plot:', error);
    plotDiv.innerHTML = `<p>Error loading or plotting Torus data: ${error.message}. Check console and ensure CSV files are generated and paths are correct.</p>`;
});

/**********************
 * HILBERT CSV LOADING
 **********************/
d3.csv('../../build/bin/locs_hilbert.csv')
  .then(function(rows) {
    function unpack(rows, key) {
      return rows.map(function(row) {
        // Ensure values are treated as numbers
        return +row[key];
      });
    }

    // Check if there's z
    var is3D = rows[0] && rows[0].hasOwnProperty('z');

    // Extract x, y, z
    var x = unpack(rows, 'x');
    var y = unpack(rows, 'y');
    var z = is3D ? unpack(rows, 'z') : null;

    // Build a color array based on index
    var colorData = [];
    for (var i = 0; i < x.length; i++) {
      colorData.push(i);
    }

    // Create the trace
    var hilbertTrace = {
      x: x,
      y: y,
      mode: 'lines+markers',
      marker: {
        size: MARKER_SIZE, // Use constant
        opacity: OPACITY, // Use constant
        color: colorData, // color each marker by index
        colorscale: COLORSCALE, // Use constant
        cmin: 0,
        cmax: x.length - 1,
        showscale: false // Hide marker color scale
      },
      line: {
        color: colorData, // color line segments by index
        colorscale: COLORSCALE, // Use constant
        cmin: 0,
        cmax: x.length - 1,
        width: LINE_WIDTH, // Use constant
        showscale: true // Show line color scale
      },
      type: is3D ? 'scatter3d' : 'scatter'
    };

    if (is3D) {
      hilbertTrace.z = z;
    }

    // Layout: 2D or 3D
    var hilbertLayout = {
      title: 'Hilbert Curve Path', // Add a title
      margin: { l: 0, r: 0, b: 0, t: 40 }, // Add top margin for title
      scene: is3D ? { // Add scene settings for 3D
          xaxis: { title: 'X', showgrid: true, zeroline: true },
          yaxis: { title: 'Y', showgrid: true, zeroline: true },
          zaxis: { title: 'Z', showgrid: true, zeroline: true },
          aspectratio: { x: 1, y: 1, z: 1 }, // Make aspect ratio cubic for 3D
          camera: { eye: {x: 1.5, y: 1.5, z: 1.5} }
      } : { // Add axis settings for 2D
          xaxis: { title: 'X', showgrid: true, zeroline: true, scaleanchor: "y", scaleratio: 1 },
          yaxis: { title: 'Y', showgrid: true, zeroline: true }
      }
    };

    Plotly.newPlot('hilbert', [hilbertTrace], hilbertLayout);
  })
  .catch(function(error) {
    console.error('Error loading or plotting Hilbert CSV file:', error);
    // Use the correct div id 'hilbert'
    document.getElementById('hilbert').innerHTML = `<p>Error loading or plotting Hilbert data: ${error.message}. Check console.</p>`;
  });