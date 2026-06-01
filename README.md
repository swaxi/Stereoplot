# qgis-stereonet v1.0.0
 WAXI QFIELD Fork of steronet plugin

# Source and Development

Original stereonet functionality developed by:

- Joe Kington (mplstereonet)
- Daniel Childs (QGIS stereonet plugin)

Major redevelopment, extension and maintenance:

- Julien Perret (Centre for Exploration Targeting, University of Western Australia)
- Mark Jessell (Centre for Exploration Targeting, University of Western Australia)



## Version 1.0.0

Version 1.0.0 represents the first major release of the WAXI/CET redevelopment of the stereonet plugin, extending the original plotting capabilities with:

- Automatic structural data recognition
- Interactive attribute classification
- Native QGIS filtering
- Dynamic contouring and girdle fitting
- Kinematic visualisation
- Rose-diagram analysis
- Persistent style templates
- Interactive category management


## Install   
Download zip file from github, install into QGIS using plugin manager   

## Usage
 1- Select a layer that has structural info in QGIS   
 2- Select the points you want to plot with one of the Select Tools (**NOT** the Identification Tool)   
 3- You can use the built in settings icon to the right of the steroenet icon OR via the WAXI QFIELD Plugin (https://github.com/swaxi/WAXI_QF) to control display behaviour   
 4- Click on WAXI Stereonet icon    ![plugin_icon](icon.png)  
    
- Planar structures can be displayed as poles or great circles   
- Linear structures are displayed as poles or rose diagrams, but if a planar feature is assocated with a linear feature, that planar feature will optionally be displayed as a great circle in a stereoplot   

## Interactive Stereonet Selection
After plotting a stereonet, you can select poles directly in the plot window and have those features highlighted in the QGIS map layer.

Note: Selection only works in poles mode (i.e., when showGtCircles is false in your config, which is the default). Great-circle and rose-diagram plots do not support selection.

#### Lasso Selection
Click and drag to draw a freehand polygon around any number of poles.

- In the stereonet window, left-click and drag to draw a lasso shape
- Release the mouse button to complete the selection
- All poles inside the lasso are highlighted with red open circles
- The corresponding features are immediately selected in the QGIS map layer
- Shift select adds points to the selection   

#### Single-Point Selection
Click near any individual pole to select it.

- Left-click close to a pole (without dragging)
- The nearest pole within the click tolerance is highlighted
- The corresponding feature is selected in the QGIS map layer
#### Clearing the Selection
Press Escape in the stereonet window to clear all selected poles and remove the selection from the QGIS map layer.

#### Tips
- The lasso and click selection replace the current selection each time — they do not add to an existing selection.
- If a pole is plotted but the click doesn't register, try clicking a little closer to the centre of the point — the tolerance is tuned to the stereonet's projection coordinate space.
- The stereonet window must remain open and in focus for keyboard shortcuts (e.g., Escape) to work.
- If best fit great circle is toggled on, the quality of the best ft is supplied as:
- **K** uses the Woodcock (1977) formula ln(e1/e2) / ln(e2/e3) where e1 ≥ e2 ≥ e3 are the covariance eigenvalues. The ±0.1 window around 1.0 labels the grey zone as "transitional" rather than forcing a binary call.
- **C** = ln(e1/e3) — it's 0 for a perfectly isotropic (random) cloud and grows without bound as the fabric tightens; there's no upper cap, so "larger = stronger" is the correct framing.
   
## Field Names

The following field names are currently recognised, you can go in the file _ _init.py__ from around line 159 and add your own if you like:

- Strike field names = ['Strike_RHR', 'Strike', 'strike']
- Dip Direction field names = ['Dip_Direction', 'Dip_Dir', 'DipDirection', 'dip_direction']
- Dip field names = ['Dip', 'dip']
- Azimuth field names = ['Azimuth', 'azimuth', 'Bearing', 'bearing']
- Plunge field names = ['Plunge', 'plunge']
- Strike of plane for lineations field names = ['Strike_ref', 'Strike_Ref', 'strike_ref']
- Dip of plane for lineations field names = ['Dip_ref', 'Dip_Ref', 'dip_ref']
- Kinematics field names = ['Kinematics', 'kinematics']
- Pitch field names = ['Pitch_RHR', 'Pitch_rhr', 'Pitch_Rhr', 'Pitch', 'pitch_rhr', 'RHR_pitch', 'rhr_pitch', 'pitch']
   
## Roadmap:




## New Features and Improvements (v1.0.0)

### Automatic Structural Data Detection
- Automatic recognition of:
  - Planes Only
  - Lineations Only
  - Lineations with Bearing Planes
- Detection is based on recognised field names rather than layer names.
- Supports both traditional planar datasets and combined planar-linear datasets.

### Bearing Plane Support
- Added support for:
  - Strike_ref + Dip_ref
  - DipDir_ref + Dip_ref
- Bearing planes are automatically displayed for combined datasets.
- Fixed initial plotting issue requiring manual checkbox toggling.

### Lineation Contouring
- Added density contouring for lineation datasets using mplstereonet line-density calculations.
- Contours now work for:
  - Lineations Only
  - Lineations with Bearing Planes

### Improved Settings Behaviour
- Automatic synchronisation between Data to Plot mode and Lineation-bearing Planes option.
- Consistent project-level persistence through stereonet.json.
- Improved fallback behaviour using QSettings.

### Kinematic Arrow Visualisation
- Added hangingwall displacement arrow plotting.
- Supports recognised kinematics fields and common naming variants.
- Supports recognised kinematic classes:
  - Sinistral
  - Dextral
  - Normal
  - Reverse / Thrust
- Validation checks ensure:
  - a valid kinematics field exists,
  - recognised kinematic values are present,
  - bearing plane information is available.

### Hangingwall Displacement Arrow Options
- User-selectable arrow construction position:
  - Plane pole
  - Lineation
- Arrows scale naturally with the stereonet during figure resizing.
- Arrow length calibrated for improved readability.
- Arrows always honour the selected kinematic sense.

### Additional Recognised Fields
- Expanded support for:
  - Trend
  - DipDir
  - DipDir_ref
  - Multiple kinematics-field aliases.
