# SMF and MMF Fiber Simulation

This repository contains MATLAB (toolbox-free) simulations of:

- Single-Mode Fiber (SMF)
- Multi-Mode Fiber (MMF)

The simulations are based on:
- Weakly guiding approximation
- LP modes
- Gaussian approximation for SMF
- Bessel-function-based modes for MMF

## Requirements
- MATLAB (base installation, no paid toolboxes)

# Project Objectives: SMF and MMF Simulation

## Step 1: Mode shapes only
- [ ] Visualize SMF Gaussian mode and MMF LP modes (Bessel + cos(l*phi))  
- [ ] Plot intensity as |E|^2  
- [ ] Purpose: understand mode shapes and superposition  

## Step 2: Include Only Physically Allowed Modes
- [ ] Compute fiber V-number: V = (2 * pi * a / lambda) * NA  
- [ ] Include only modes with u_lm <= V  
- [ ] Purpose: reflect which modes the fiber can physically support  

## Step 3: Normalize and Weight Modes
- [ ] Assign amplitude coefficients c_lm based on input beam (centered, offset, tilted)  
- [ ] Normalize each mode to unit maximum or unit power  
- [ ] Purpose: realistic mode excitation  
- [ ] Effect: the mode superposition resembles real excitation  

## Step 4: Include Propagation Phase (beta)
- [ ] Each mode propagates along z: E_lm(r,phi,z) = E_lm(r,phi) * exp(i * beta_lm * z)  
- [ ] Compute beta_lm approximately as k0 * n_eff  
- [ ] Sum modes at each z to get interference pattern  
- [ ] Purpose: visualize changing speckle along fiber length  

## Step 5: Include Cladding / Evanescent Field
- [ ] Replace the mask (R <= a) with proper cladding solution using modified Bessel functions  
- [ ] Purpose: more physically accurate mode confinement  

## Step 6: Modal Coupling (Advanced)
- [ ] Introduce mode coupling due to imperfections or bending:  
  - dc_lm/dz = sum over l',m' of (kappa_lm,l'm' * c_l'm')  
- [ ] Purpose: realistic power redistribution between modes  

## Step 7: Time-Dependent Simulation (Optional)
- [ ] Add temporal or wavelength dependence  
- [ ] Simulate pulse propagation and modal dispersion  
- [ ] Purpose: observe pulse broadening and speckle evolution  

## Step 8: Visualization Improvements
- [ ] 2D intensity: imagesc(abs(E).^2)  
- [ ] 3D surface: surf(X,Y,abs(E).^2)  
- [ ] Animate propagation along z for speckle patterns  

