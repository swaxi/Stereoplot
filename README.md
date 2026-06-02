# Stereoplot Plugin v1.0.0

## Overview

The Stereoplot Plugin is an interactive stereographic projection plotting and analysis tool for QGIS.

The plugin supports:

- Stereonet plotting of planar and linear structures
- Rose diagram generation
- Density contouring (Modified Kamb method)
- Best-fit girdle analysis
- Interactive stereonet selection linked to QGIS
- Bearing-plane visualisation
- Kinematic arrow plotting (hangingwall displacement for S-L fabrics)
- Attribute-based classification
- Attribute filtering using native QGIS expressions
- Category-specific styling
- Style template management
- Publication-ready figure export

Version 1.0.0 represents the first major release of the WAXI/CET redevelopment of the original stereonet plugin.

---

# Installation

1. Download the latest ZIP package from GitHub.
2. Open **QGIS → Plugins → Manage and Install Plugins**.
3. Select **Install from ZIP**.
4. Browse to the downloaded ZIP file.
5. Install and enable the plugin.

---

# Quick Start

## Basic Workflow

### Step 1 – Select a Layer

Select a point layer containing structural measurements.

The plugin automatically detects whether the layer contains:

- Planar measurements
- Linear measurements
- Lineations with associated bearing planes

### Step 2 – Select Features

Select the features you wish to plot using any standard QGIS selection tool.

**Important:** The Identify Tool does not create a selection and cannot be used for plotting.

### Step 3 – Configure Settings (Optional)

Open the settings dialog using:

- the Stereoplot Settings button, or
- the WAXI QField Plugin.

### Step 4 – Generate Plot

Click the Stereoplot button.

The plugin will automatically determine the appropriate plotting mode.

---

# Supported Data Types

## Planes Only

Requires:

- Strike + Dip

or

- Dip Direction + Dip

Displayed as:

- poles to planes
- great circles

---

## Lineations Only

Requires:

- Trend/Azimuth/Bearing
- Plunge

Displayed as:

- lineations
- rose diagrams
- density contours

---

## Lineations with Bearing Planes

Requires:

- Trend/Azimuth/Bearing
- Plunge

and either:

- Strike_ref + Dip_ref

or

- DipDir_ref + Dip_ref

Displayed as:

- lineations
- associated bearing planes
- kinematic arrows
- density contours

---

# Plot Types

## Stereonet

Supports:

- Poles to planes
- Great circles
- Lineations
- Bearing planes
- Kinematic arrows
- Density contours
- Best-fit girdles

---

## Rose Diagram

Rose diagrams may be generated from:

- planar datasets
- linear datasets

Features:

- automatic orientation extraction
- dynamic resizing
- publication-ready output

---

# Density Contours

Density contouring uses the Modified Kamb method implemented through mplstereonet.

Supported for:

- poles to planes
- lineations
- lineations with bearing planes

Features:

- continuous colour scale
- σ-level contour labels
- dynamic updates during classification and filtering
- colour-bar export

---

# Best-Fit Girdle Analysis

Best-fit girdles can be calculated from the currently visible dataset.

The girdle automatically updates when:

- categories are hidden or shown
- filters are modified
- classifications are changed

### Fabric Statistics

#### Woodcock K-value

Computed as:

K = ln(e1/e2) / ln(e2/e3)

where:

e1 ≥ e2 ≥ e3

are covariance eigenvalues.

#### Fabric Strength (C)

Computed as:

C = ln(e1/e3)

Higher values indicate stronger fabrics.

---

# Kinematic Arrow Visualisation

## Purpose

Kinematic arrows display the inferred hangingwall displacement direction.

---

## Supported Kinematic Classes

Recognised values include:

### Sinistral

- Sinistral
- Sinistral-slip
- Left-lateral
- Sin

### Dextral

- Dextral
- Dextral-slip
- Right-lateral
- Dex

### Normal

- Normal
- Normal-slip
- Extensional

### Reverse

- Reverse
- Reverse-slip
- Thrust
- Compressional

---

## Kinematics Fields

Common recognised field names include:

- Kinematics
- Kinematic
- Kin
- Movement
- SlipSense
- Slip_Sense
- ShearSense
- SenseOfMovement

Field matching is case-insensitive.

---

## Hangingwall Displacement Arrow Position

Users may choose whether arrows are constructed from:

- Plane Pole
- Lineation

This option is available from the Settings dialog.

---

# Attribute Classification

## Purpose

Classification allows structural data to be grouped according to attribute values.

Examples:

| Field | Categories |
|---------|---------|
| Generation | 0, 1, 2 |
| Kinematics | Sinistral, Dextral |
| Lithology | Mafic, Felsic |

---

## Enabling Classification

1. Enable **Classification**.
2. Select a classification field.
3. Click **Update Settings**.

The plugin automatically creates categories from unique attribute values.

---

# Category Management

For each category users can:

- Show/Hide
- Select All
- Hide All
- Invert Selection

Updates occur dynamically without recreating the stereonet.

---

# Category Styling

Each category may have independent:

- Symbol shape
- Symbol colour
- Symbol size
- Line colour
- Line width
- Transparency
- Arrow colour

Changes are applied immediately.

---

# Style Management

Style templates can be:

- Saved
- Loaded
- Reset
- Deleted

Templates are stored in:

```text
stereonet_styles.json
```

within the project configuration directory.

---

# Attribute Filtering

## Purpose

Filter records before plotting.

## Filter Builder

The plugin uses the native QGIS expression system.

Examples:

```sql
"Generation" = '1'
```

```sql
"Kinematics" IN ('Sinistral-slip','Reverse-slip')
```

```sql
"Generation" = '1'
AND "Kinematics" IS NOT NULL
```

Supported:

- AND
- OR
- NOT
- Numerical comparisons
- Text comparisons
- NULL checks
- Nested expressions

---

# Interactive Selection

Selections made within the stereonet are synchronised with QGIS.

Supported in:

- Pole plots

Not supported in:

- Great-circle plots
- Rose diagrams

## Lasso Selection

Click and drag around points.

Selected features are:

- highlighted on the stereonet
- selected in QGIS

## Single Point Selection

Click near a plotted point.

The nearest feature is selected.

## Clear Selection

Press:

```text
Esc
```

to clear selections.

---

# Exporting Figures

Figures can be exported directly from the plot window.

Exported figures include:

- stereonet
- contours
- colour scales
- best-fit girdles
- classification legend
- rose diagrams
- kinematic arrows

Suitable for publication and reporting.

---

# Recognised Structural Fields

## Planes

### Strike

- Strike_RHR
- Strike
- strike

### Dip Direction

- Dip_Direction
- Dip_Dir
- DipDirection
- DipDir

### Dip

- Dip
- dip

## Lineations

### Trend / Azimuth

- Azimuth
- Bearing
- Trend

### Plunge

- Plunge

## Bearing Planes

### Strike

- Strike_ref
- Strike_Ref

### Dip Direction

- DipDir_ref
- DipDirection_ref

### Dip

- Dip_ref
- Dip_Ref

---

# Credits

## Original Development

- Joe Kington (mplstereonet)
- Daniel Childs (QGIS Stereonet)

## WAXI/CET Redevelopment

- Julien Perret — Centre for Exploration Targeting, University of Western Australia
- Mark Jessell — Centre for Exploration Targeting, University of Western Australia

---

# Version

Current version: **1.0.0**
