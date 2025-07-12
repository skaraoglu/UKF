# Unscented Kalman Filter (UKF) for Parameter Estimation in Nonlinear Dynamical Systems

This R package implements the Unscented Kalman Filter (UKF) for parameter estimation and state tracking in nonlinear dynamical systems, with a focus on applications such as fMRI time series analysis and coupled oscillator models.

---

## Overview

The Unscented Kalman Filter (UKF) is a recursive Bayesian filter designed for nonlinear systems. Unlike the Extended Kalman Filter (EKF), which linearizes the system using Jacobians, the UKF uses a deterministic sampling technique (the unscented transform) to more accurately capture the mean and covariance of the system state under nonlinear transformations.

This package provides tools for:
- State and parameter estimation in nonlinear ODE models
- Iterative optimization of model parameters using UKF
- Visualization and smoothing of time series data

---

## Background

The Unscented Kalman Filter (UKF) is a recursive Bayesian filter designed for nonlinear dynamical systems. It is particularly well-suited for systems where the state evolution and/or observation models are nonlinear, and it avoids the need for explicit Jacobian calculations required by the Extended Kalman Filter (EKF).

*Nonlinear State-Space Model*

The general discrete-time nonlinear state-space model is:

$$
\begin{aligned}
\mathbf{x}_{k+1} &= f(\mathbf{x}_k) + \mathbf{w}_k \\
\mathbf{y}_k &= h(\mathbf{x}_k) + \mathbf{v}_k
\end{aligned}
$$

where:

- $\mathbf{x}_k$ is the state vector at time $k$
- $\mathbf{y}_k$ is the observation vector
- $f(\cdot)$ is the (possibly nonlinear) state transition function
- $h(\cdot)$ is the (possibly nonlinear) observation function
- $\mathbf{w}_k \sim \mathcal{N}(0, Q)$ is the process noise
- $\mathbf{v}_k \sim \mathcal{N}(0, R)$ is the measurement noise

Note: In this package, the state vector is augmented to include both the system states and unknown model parameters, allowing for simultaneous parameter and state estimation.

### UKF Process

The UKF operates in two main phases at each time step: prediction and update. Here is a detailed breakdown of the code:

1. Sigma Point Generation
   - Given the current estimate of the augmented state mean $\hat{\mathbf{x}}{k-1}$ and covariance $P{k-1}$, generate $2L$ sigma points (where $L$ is the dimension of the augmented state).
   - The sigma points are chosen to capture the mean and covariance of the state distribution. 
2. Propagation (Prediction Step)
   - Each sigma point is propagated through the nonlinear system dynamics.

   - In this implementation, this is done using a Runge-Kutta 4th order (RK4) integrator for the ODE model. This provides a numerically stable and accurate way to propagate the state and parameter estimates forward in time, even for stiff or highly nonlinear systems.

   - The predicted mean and covariance are computed from the propagated sigma points.

   - The use of the Runge-Kutta 4th order (RK4) method for propagating sigma points is a key feature of this implementation. RK4 is a widely used, robust numerical integrator for ODEs. It provides a good balance between accuracy and computational efficiency, making it suitable for the nonlinear ODEs encountered in biological and physical systems.

3. Measurement Prediction
   - The propagated sigma points are mapped through the observation function $h(\cdot)$ to predict the expected measurement for each sigma point.
   - The predicted measurement mean and covariance are computed.
4. Update (Correction Step)
   - The cross-covariance between the state and measurement is computed.
   - The Kalman gain is calculated, which determines how much the state estimate should be corrected based on the new measurement.
   - The state mean and covariance are updated using the actual measurement.
5. Process and Measurement Noise
   - Process noise ($Q$): Added to the parameter block of the covariance matrix at each prediction step, allowing for uncertainty in parameter evolution.
   - Measurement noise ($R$): Used in the update step to account for observation uncertainty.

---

## UKF Methods: Streamlined vs. Standard

This package provides two UKF one-step methods for parameter and state estimation:

- **Streamlined UKF (`UKF_dT`)**:  
  - Uses $2N_x$ sigma points (no central mean point or scaling parameter).
  - All sigma points are equally weighted in mean and covariance calculations.
  - Covariance update uses $P_{xx} = P_{xx} - K P_{xy}^T$ (empirically robust for this application).
  - Process noise $Q$ is only added to the parameter block of the covariance matrix.
  - Designed for computational efficiency and robust optimization in moderate state dimensions.

- **Standard UKF (`UKF_dT_std`)**:  
  - Implements the original UKF formulation with $2L+1$ sigma points and scaling parameter $\lambda$.
  - Uses UKF weights ($W_m^i$, $W_c^i$) for mean and covariance.
  - Covariance update follows $P_{xx} = P_{xx} - K P_{yy} K^T$.
  - Robust positive-definite covariance handling using the `Matrix` package.
  - Suitable for users who want strict adherence to the original UKF algorithm.

Both methods propagate sigma points using a Runge-Kutta 4th order (RK4) integrator for nonlinear ODE models.  
The observation model is assumed linear (direct selection of state variables).

**Select the method via the `method` parameter in `UKF_blend` and `iterative_param_optim` (`"streamlined"` or `"standard"`).**

---

## Additional Parameters

- `param_tol`: Tolerance for parameter convergence (default: `1e-4`)
- `MAXSTEPS`: Maximum number of UKF iterations (default: `30`)
- `forcePositive`: If `TRUE`, enforces positivity for parameters
- `seeded`: If `TRUE`, sets random seed for reproducibility
- `method`: `"streamlined"` or `"standard"` UKF step

---

## Installation

```r
library(devtools)
install_github("insilico/UKF")
library(UKF)
```

---

## Usage Example

See `inst/model_voxels_2param.R` for a full example.  
Here is a minimal workflow:

```r
# Prepare your data: first column is time, others are observed variables
data <- read.csv("data/your_timeseries.csv")

# Define your ODE model as a function
my_ode_model <- function(t, y, p) {
  # y: state vector, p: parameter vector
  # return dydt as a vector
}

# Set initial parameter guess and UKF settings
param_guess <- c(1, 1) # example for two parameters
N_p <- 2
N_y <- ncol(data) - 1
dt <- 0.1
dT <- 1.0

# Run UKF parameter estimation (choose method: "streamlined" or "standard")
result <- iterative_param_optim(
  param_guess, t_dummy = 0, ts_data = data, ode_model = my_ode_model,
  N_p = N_p, N_y = N_y, dt = dt, dT = dT,
  param_tol = 1e-4, MAXSTEPS = 100,
  R_scale = 0.3, Q_scale = 0.01,
  trace = TRUE, forcePositive = TRUE,
  seeded = TRUE,
  method = "standard" # or "streamlined"
)

# Access estimated parameters and diagnostics
result$param_est
result$chisq
plot_ukf_chi_square_loss(result)
```

---

## Tuning Q and R

- `Q_scale`: Controls process noise (model/parameter uncertainty). Higher values allow parameters to adapt more quickly but may introduce noise.
- `R_scale`: Controls measurement noise. Higher values make the filter trust the model more than the data.

**Tip:** Start with `Q_scale = 0.01` and `R_scale = 0.3`, then adjust based on the stability and fit of your results.

---

## Visualization

The package includes plotting functions for:
- Raw and smoothed time series (`plot_voxels_and_smooth_ggplot`)
- UKF state and parameter estimates (`plot_ukf_and_smoothed_ggplot`)
- Chi-square loss over iterations (`plot_ukf_chi_square_loss`)

---

## Dependencies

Some functions require:
```r
install.packages(c('KernSmooth', 'ggplot2', 'pracma', 'Matrix'))
```

---

## Academic Use

This package is developed for research in nonlinear time series analysis, parameter estimation, and model validation.  
If you use this package in your research, please cite appropriately and contact the maintainer for collaboration or questions.

---

## Contact

For questions, bug reports, or collaboration, contact:  
[brett.mckinney@gmail.com](mailto:brett.mckinney@gmail.com)