<h1 align="center">Closed-Form Construction of Voronoi Diagrams with Star-Shaped Metrics</h1>

<p align="center">
  <img src="images/teaser.svg" alt="Teaser image showing Voronoi diagrams with star-shaped metrics on surfaces." style="width:100%;height:auto;display:block;margin:0 auto;max-width:1600px;" />
</p>

<div align="center">
<a href="mailto:haoyang.zhou@ai.ethz.ch">Haoyang Zhou</a>, <a href="mailto:logan.numerow@inf.ethz.ch">Logan Numerow</a>, <a href="mailto:scoros@inf.ethz.ch">Stelian Coros</a>, <a href="mailto:bthomasz@ethz.ch">Bernhard Thomaszewski</a>
<br/>
ETH Zürich, Switzerland
</div>

----

This repository contains the implementation for the paper "Closed-Form Construction of Voronoi Diagrams with Star-Shaped Metrics" published in ACM Transactions on Graphics and presented in SIGGRAPH Asia 2025. The code provides algorithms, utilities, and example programs to construct Voronoi diagrams induced by star-shaped metrics both in the 2D plane and on surfaces.

<div align="center">
<a href="https://doi.org/10.1145/3763296">📑 Paper</a>
<a href="https://ssmvd.haoyang.io/">🕹️ Online Demo</a>
<a href="https://www.icloud.com/keynote/036lErjBS7y9mGKStJPA3NgGw#SA25-Slides-Final">🎞️ Presentation Slides</a>
</div>

## Getting Started

### Contents
- `include/` — public headers for the library
- `examples/` — example programs and demos (see individual subfolders)
  - `Surface_Voronoi_diagram_with_star_metrics/` 
    - `random_sites.cpp` — build and visualize a diagram on surfaces with random sites
  - `Voronoi_diagram_with_star_metrics_2/`
    - `random_sites.cpp` — build and visualize a diagram with random sites
    - `grid_sites.cpp` — build and visualize a regular diagram with grid sites
- `tests/` — test helpers and small test programs
- `CMakeLists.txt` — top-level CMake project

### Build

#### Prerequisites
- **C++ toolchain**: A standards-compliant C++ compiler (C++20 or newer).
- **CMake**: version 3.16+.
- **Libraries**: CGAL and its [dependencies](https://doc.cgal.org/latest/Manual/thirdparty.html) (Boost and Qt6 for visualization in examples).

Refer to [CGAL installation instructions](https://doc.cgal.org/latest/Manual/general_intro.html) for help setting up CGAL and its dependencies. This project also provides options to fetch CGAL automatically with the CMake option `SSM_INCLUDE_CGAL`.

#### Build

To build the project with examples, clone the repository and build with CMake:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -B build -S .
cmake --build build
```

The build examples will be available in the `examples/` subfolders of the build output directory, with the same structure as in the source `examples/` folder.

#### Using the library from your code

Add this repository as a subdirectory in your CMake project or fetch it with `FetchContent`. Then link to the target `ssm` in your executable or library target. Please refer to the source code in the `examples/` folder for usage patterns.

## Citation

If you use this code in published work, please consider citing the paper:

Haoyang Zhou, Logan Numerow, Stelian Coros, and Bernhard Thomaszewski. 2025. Closed-Form Construction of Voronoi Diagrams with Star-Shaped Metrics. *ACM Trans. Graph.* 44, 6, Article 235 (December 2025), 13 pages. https://doi.org/10.1145/3763296

<details>
<summary>BibTex</summary>
  
```bib
@article{10.1145/3763296,
author = {Zhou, Haoyang and Numerow, Logan and Coros, Stelian and Thomaszewski, Bernhard},
title = {Closed-Form Construction of Voronoi Diagrams with Star-Shaped Metrics},
year = {2025},
issue_date = {December 2025},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {44},
number = {6},
issn = {0730-0301},
url = {https://doi.org/10.1145/3763296},
doi = {10.1145/3763296},
abstract = {Cellular patterns, from planar ornaments to architectural surfaces and mechanical metamaterials, blend aesthetics with functionality. Homogeneous patterns like isohedral tilings offer simplicity and symmetry but lack flexibility, particularly for heterogeneous designs. They cannot smoothly interpolate between tilings or adapt to double-curved surfaces without distortion. Voronoi diagrams provide a more adaptable patterning solution. They can be generalized to star-shaped metrics, enabling diverse cell shapes and continuous grading by interpolating metric parameters. Mart\'{\i}nez et al. [2019] explored this idea in 2D using a rasterization-based algorithm to create compelling patterns. However, this discrete approach precludes gradient-based optimization, limiting control over pattern quality. We introduce a novel, closed-form, fully differentiable formulation for Voronoi diagrams with piecewise linear star-shaped metrics, enabling optimization of site positions and metric parameters to meet aesthetic and functional goals. It naturally extends to arbitrary dimensions, including curved 3D surfaces. For improved on-surface patterning, we propose a per-sector parameterization of star-shaped metrics, ensuring uniform cell shapes in non-regular neighborhoods. We demonstrate our approach by generating diverse patterns, from homogeneous to continuously graded designs, with applications in decorative surfaces and metamaterials.},
journal = {ACM Trans. Graph.},
month = dec,
articleno = {236},
numpages = {13},
keywords = {voronoi diagram, tiling, pattern, differentiable simulation, metamaterial}
}
```
</details>

## Contact & Contributions
- **Bug reports / issues**: Please open an issue on this repository with reproduction steps and any build logs.
- **Contributions**: Contributions are welcome. Please open a pull request and follow the existing code style and CMake structure.
- **Authors / contact**: For questions about the algorithm or paper details, contact the paper authors (listed in the citation).
