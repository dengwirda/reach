## `REACH: Stream-aligned unstructured mesh generation`

Generate variable resolution unstructured meshes for large-scale geophysical domains with vector-based hydrological features 'burned-in' to the mesh topology.

<p align="middle">
  <img src = "../master/img/reach_msh.png" hspace="0.25%">
</p>

### `Quickstart`

    Download the HydroSHEDS-v1 shapefile (HydroRIVERS) to reach/dat
    Build the cython back-end: python3 setup.py build_ext --inplace
    Run a simple example case: python3 example.py

The HydroSHEDS-v1 vector database is available [here](https://www.hydrosheds.org/products/hydrorivers).

