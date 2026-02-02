# MATLAB code

## Dependencies
- Manopt (Riemannian optimization in MATLAB): https://www.manopt.org
- Optional for global certification:
  - GloptiPoly 3: https://homepages.laas.fr/henrion/software/gloptipoly3/
  - An SDP solver supported by GloptiPoly such as SeDuMi.

## Quick start
1. Add Manopt to your MATLAB path.
2. Run the continuous-motion simulation and product-manifold solver:

```matlab
addpath(genpath('path/to/manopt'));
addpath(genpath(pwd));
demo_manopt_product_manifold
```

3. Optional global SOS windowed solve:

```matlab
demo_global_sos_gloptipoly
```

The simulation generates continuous (dense) magnetometer readings by integrating a smooth angular velocity and sampling at a fixed rate.

## Real data: UST Nautilus TXT logs

This package also includes loaders for the UST Nautilus log format:

Each line contains 10 numeric fields:
`seq gx gy gz ax ay az mx my mz`.

The magnetometer samples are the last three fields.

Run the real-data experiment (requires Manopt; GloptiPoly optional):

```matlab
addpath(genpath('path/to/manopt'));
addpath(genpath(pwd));
run('realdata_nautilus/run_nautilus_experiment.m');
```
