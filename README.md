# MATLAB codes for a projection method for eigenvalue problems of linear nonsquare matrix pencils

by Keiichi Morikuni, Ph.D.

The program is designed to compute the eigenvalues in a circle region and the corresponding eigenvectors for a matrix pencil.

## Usage

To generate a tall test matrix pencil, execute the following:
```
>> scrpt4tall
```
or to generate a wide test matrix pencil, execute the following:
```
>> scrpt4wide
```

To run the main function, execute the following:
```
>> scrpt4driver
```
The elapsed CPU time, relative residual norm, and the relative error of the computed eigenvalues are output, as well as the plots of the exact eigenvalues and computed eigenvalues in the complex plane.

<img src="https://user-images.githubusercontent.com/15831262/137937710-319a4a88-c14f-4702-a52e-418b2c5567f4.jpg" width="520pt">

### Reference
Keiichi Morikuni, Projection method for eigenvalue problems of linear nonsquare matrix pencils,
SIAM Journal on Matrix Analysis and Applications, Volume 42, Number 3,
pp. 1381-1400, September 20, 2021. DOI: [10.1137/20M1377886](https://doi.org/10.1137/20M1377886)

If you use the codes in research for publication, please cite this paper.
```bibtex
@Article{Morikuni2021SIMAX,
  author    = {Keiichi Morikuni},
  doi       = {10.1137/20m1377886},
  journal   = {SIAM Journal on Matrix Analysis and Applications},
  number    = {3},
  pages     = {1381--1400},
  title     = {Projection method for eigenvalue problems of linear nonsquare matrix pencils},
  volume    = {42},
  year      = {2021},
  publisher = {Society for Industrial {\&} Applied Mathematics ({SIAM})},
}
```
