# 3D-CFD-with-Weissinger-method

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=lucapiombo/3D-CFD-with-Weissinger-method)

This project aims to perform a comprehensive analysis of a three-dimensional aerodynamic problem by examining a finite wing through the numerical implementation of Weissinger’s method.

The primary objective is to design a finite wing that minimizes induced drag while maintaining a fixed angle of attack and achieving a specified range of lift coefficient values (cL). Additionally, the effect of a tail is considered to reduce the wing's pitching moment, with analyses performed both with and without ground effects. The methodology has been implemented in MATLAB, and its validation has been carried out using XFLR5.

## Project Structure

The repository contains several files along with three main scripts, each corresponding to a different configuration:
- **Entire Wing Analysis:** Computes results for the entire wing.
- **Tandem Configuration:** Analyzes a tandem wing configuration.
- **Tandem with Ground Effect:** Evaluates a tandem configuration while considering ground effects.

## Getting Started

1. Clone this repository.
2. Open MATLAB and navigate to the project directory.
3. Run the appropriate main script to generate and view the results for your desired configuration.

---

Feel free to modify or extend the code for further analysis or research.
