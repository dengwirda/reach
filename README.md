## `REACH: Stream-aligned unstructured mesh generation`

Generate filtered stream networks and variable resolution unstructured meshes for large-scale domains with embedded vector-based hydrological features 'burned-in' to the mesh topology.

<p align="middle">
  <img src = "../main/img/reach_msh.png" hspace="0.25%">
</p>

### `Quickstart`

    Download the HydroSHEDS-v1 shapefile (HydroRIVERS) to reach/dat
    Install dependencies: pip3 install -r requirements.txt
    Build the cython back-end: python3 setup.py build_ext --inplace
    Run example cases: python3 example.py

The required HydroSHEDS-v1 vector database is available [here](https://www.hydrosheds.org/products/hydrorivers).

