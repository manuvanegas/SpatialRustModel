[![DOI](https://zenodo.org/badge/677233131.svg)](https://zenodo.org/badge/latestdoi/677233131)

# SpatialRust Model

The SpatialRust model is a computational tool designed to simulate the dynamics of coffee leaf rust (CLR) epidemics in agroforestry systems. This model integrates coffee plant physiology, shade tree dynamics, and CLR progression to provide insights into the complex interactions affecting CLR epidemics and coffee production.

The model was developed as the center of my doctoral dissertation in Biological Design. It does not produce quantitative predictions, but it is a versatile tool to understand trends and establish causal links between the components of this complex system. I will add a link to the dissertation document once it is made available.

## Key Features

- **Integrated Approach:** SpatialRust integrates coffee plant growth, shade tree dynamics, and CLR epidemic processes to provide a comprehensive view of how these factors interact.

- **Spatially Explicit:** The model accounts for spatial variability in shade distribution, allowing for a more realistic representation of the microclimate within a coffee agroforestry system.

- **Management Strategies:** The model allows users to explore different management strategies for CLR control, considering factors such as shading patterns and farm management practices.

- **Multi-Dimensional Insights:** By simulating the interactions among coffee plants, shade trees, and CLR, the model provides insights into how different combinations of factors can impact coffee production and disease progression.

## Usage

1. Install the required dependencies. Move to this project's directory and run:
```bash
$ julia scripts/install.jl
```
2. Run the model using the provided sample script.
```bash
$ julia scripts/samplerun.jl
```
4. Adjust parameters in the script to explore different scenarios and management strategies (a description of the different parameters will be added here once the dissertation document is made available).
5. Review the model's outputs (in the `results` folder).

## License

The SpatialRust model is released under the MIT license. Please review the LICENSE file for more details.

## Contact

If you have any questions, suggestions, or feedback, please [create an issue](https://github.com/manuvanegas/SpatialRustModel/issues).

## Citation

A manuscript is in preparation. Check back later or reach out through the repository's [issues](https://github.com/manuvanegas/SpatialRustModel/issues) if you plan to use the model!

