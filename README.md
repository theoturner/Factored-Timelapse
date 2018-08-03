# Factored Timelapse
### Scene lighting componentisation for advanced image editing

Time-lapse photography captures a large amount of information about the lighting of a scene over a period of time.
If this data is factored into its individual components – illumination, reflectance and shadows – the scene may be computationally re-lit for the purposes of advanced image editing operations.
This program implements the Sunkavalli et al. 2007 approach to factored timelapses.

### Outline
The algorithm involves splitting scene lighting into sunlight and skylight components:
1.  Read time-lapse video and mask sky pixels.
2.  For each colour channel, independently:
    1. Detect shadows.
    2. Extract skylight and the skylight basis curve.
    3. Extract sunlight and the sunlight basis curve.
3.  Reconstruct images and visualise results.

An in-depth description of how each component functions is detailed in the `report.pdf` file.
