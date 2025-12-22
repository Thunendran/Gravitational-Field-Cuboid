%% ============================================================
%   Example 01 — Gravitational Potential at a Single Point
%   MATLAB Version
% ============================================================

clear; clc;

%% Add utils_MATLAB folder to path
addpath('/Volumes/Dunendran/Programs/MATLAB/Cuboid_Extension/utils_MATLAB');

%% Physical constants
G   = 1.0;
rho = 1.0;

%% Cuboid center & dimensions
X0 = 0.0;   Y0 = 0.0;   Z0 = 0.0;
L  = 1.0;   B  = 1.0;   D  = 1.0;

% MATLAB uses half-sizes
Lh = L/2;   Bh = B/2;   Dh = D/2;

%% Evaluation point
x = 0.0;
y = 0.0;
z = 0.0;

%% Compute potential using the class
V = GravCuboidPotential.potential_single( ...
        x - X0, ...
        y - Y0, ...
        z - Z0, ...
        Lh, Bh, Dh, ...
        rho, G );

%% Print result
fprintf("The gravitational potential at point (%.3f, %.3f, %.3f) is %.13e m^2/s^2\n", ...
    x, y, z, V);
