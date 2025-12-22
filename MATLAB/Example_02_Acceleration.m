%% ============================================================
% Example 02 — Validate gravitational acceleration
% ============================================================

clear; clc;

%% Add utils to path
addpath('/Volumes/Dunendran/Programs/MATLAB/Cuboid_Extension/utils_MATLAB');

%% ---------------------------------------------------------
% PARAMETERS
% ---------------------------------------------------------
a = 1.0;
b = 1.0;
c = 1.0;

rho = 1.0;
G   = 1.0;

% Point of evaluation
x = 0.0;
y = 0.0;
z = 0.0;

h = 1e-6;

results = struct();


%% ============================================================
% g_x COMPARISON
% ============================================================

gx_analytic = GravCuboidAcceleration.accel_gx(x, y, z, a, b, c, rho, G);

V_plus_x  = GravCuboidPotential.potential_single(x + h, y, z, a, b, c, rho, G);
V_minus_x = GravCuboidPotential.potential_single(x - h, y, z, a, b, c, rho, G);

gx_numeric = -(V_plus_x - V_minus_x) / (2*h);

results.gx = [gx_analytic, gx_numeric, gx_analytic - gx_numeric];


%% ============================================================
% g_y COMPARISON
% ============================================================

gy_analytic = GravCuboidAcceleration.accel_gy(x, y, z, a, b, c, rho, G);

V_plus_y  = GravCuboidPotential.potential_single(x, y + h, z, a, b, c, rho, G);
V_minus_y = GravCuboidPotential.potential_single(x, y - h, z, a, b, c, rho, G);

gy_numeric = -(V_plus_y - V_minus_y) / (2*h);

results.gy = [gy_analytic, gy_numeric, gy_analytic - gy_numeric];


%% ============================================================
% g_z COMPARISON
% ============================================================

gz_analytic = GravCuboidAcceleration.accel_gz(x, y, z, a, b, c, rho, G);

V_plus_z  = GravCuboidPotential.potential_single(x, y, z + h, a, b, c, rho, G);
V_minus_z = GravCuboidPotential.potential_single(x, y, z - h, a, b, c, rho, G);

gz_numeric = -(V_plus_z - V_minus_z) / (2*h);

results.gz = [gz_analytic, gz_numeric, gz_analytic - gz_numeric];


%% ============================================================
% PRINT RESULTS
% ============================================================

fprintf("=====================================================\n");
fprintf("  Analytical vs Numerical Gravitational Acceleration\n");
fprintf("               at (0.0, 0.0, 0.0) for a Cube\n");
fprintf("=====================================================\n");

fprintf("\n--- g_x Component ---\n");
fprintf("Analytical g_x   = %.15e\n", results.gx(1));
fprintf("Numerical g_x    = %.15e\n", results.gx(2));
fprintf("Difference       = %.3e\n",   results.gx(3));

fprintf("\n--- g_y Component ---\n");
fprintf("Analytical g_y   = %.15e\n", results.gy(1));
fprintf("Numerical g_y    = %.15e\n", results.gy(2));
fprintf("Difference       = %.3e\n",   results.gy(3));

fprintf("\n--- g_z Component ---\n");
fprintf("Analytical g_z   = %.15e\n", results.gz(1));
fprintf("Numerical g_z    = %.15e\n", results.gz(2));
fprintf("Difference       = %.3e\n",   results.gz(3));

fprintf("=====================================================\n");
