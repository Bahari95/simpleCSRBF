# simpleCSRBF
accelerated minimalist library for the compactly supported radial basis function CSRBF to simulate various problems in 1d&2D

### These codes are starting points for the more general solver.(under development)

Below is an example of Poisson equation with homogeneous Dirichlet boundary condition


#### For your analysis You can find and install the following packages from the url

# For pyccel :
  
  https://github.com/pyccel/pyccel
## Install

**Standard mode**

```shell
python3 -m pip install .
```

**Development mode**

```shell
python3 -m pip install --user -e .


# ... (left) adapted mesh (right) approximate solution
![PNG](https://github.com/Bahari95/simpleCSRBF/blob/main/r_refinement_ex.png)

# ... Stencil Matrix for adapted mesh
![PNG](https://github.com/Bahari95/simpleCSRBF/blob/main/r-refinement_matrix.png)

# ... Stencil Matrix for uniform mesh
![PNG](https://github.com/Bahari95/simpleCSRBF/blob/main/uniform_matrix.png)
