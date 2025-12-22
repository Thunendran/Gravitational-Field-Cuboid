clc; clear;

% ---------------------------------------------------------
% CONSTANTS
% ---------------------------------------------------------
G = 1.0;
rho = 1.0;

X_0 = 3.0; 
Y_0 = 3.0; 
Z_0 = 3.0;

L_dim = 2.0; 
B_dim = 2.0; 
D_dim = 2.0;

epsilon = 1e-10;

% Semi-dimensions
a = L_dim/2;
b = B_dim/2;
c = D_dim/2;

% ---------------------------------------------------------
% TEST POINTS 
% ---------------------------------------------------------
test_points = {
    [3.5, 3.5, 3.5],            "Interior (centered below face ABCD)";
    [3.5, 3.5, 4.0],            "On face ABCD";
    [3.5, 3.5, 4.0 + epsilon],  "Exterior (above face)";
    [3.5, 3.5, 4.0 - epsilon],  "Interior (below face)";
    [4.0, 4.0, 4.0],            "On vertex B";
    [4.0, 4.0 + epsilon, 4.0],  "On extended edge AB";
    [4.0 + epsilon, 4.0, 4.0],  "On extended edge CB";
    [3.0, 5.0, 4.0],            "On extended face ABCD";
    [7.0, 7.0, 7.0],            "Far exterior point";
    [7.0 + epsilon, 7.0 + epsilon, 7.0 + epsilon], "Far exterior + eps"
};

numTests = size(test_points, 1);

% Storage for table printing
results = struct([]);

% ---------------------------------------------------------
% EVALUATION LOOP (analytic tensor)
% ---------------------------------------------------------
for k = 1:numTests
    p_abs = test_points{k,1};
    desc  = test_points{k,2};

    xa = p_abs(1); 
    ya = p_abs(2); 
    za = p_abs(3);

    % Convert to cuboid-centered coordinates
    xr = xa - X_0;
    yr = ya - Y_0;
    zr = za - Z_0;

    % Analytic second derivatives
    Vxx = GravCuboidTensor.Vxx(xr, yr, zr, a, b, c, rho, G);
    Vyy = GravCuboidTensor.Vyy(xr, yr, zr, a, b, c, rho, G);
    Vzz = GravCuboidTensor.Vzz(xr, yr, zr, a, b, c, rho, G);

    Lap = Vxx + Vyy + Vzz;

    results(k).AbsPoint = sprintf("(%.10f, %.10f, %.10f)", xa, ya, za);
    results(k).Description = desc;
    results(k).Vxx = Vxx;
    results(k).Vyy = Vyy;
    results(k).Vzz = Vzz;
    results(k).Laplacian = Lap;
end

% ---------------------------------------------------------
% PRINT RESULTS
% ---------------------------------------------------------
fprintf("Gravitational Tensor Laplacian Test Results:\n");
fprintf("Cube Center = (%.1f, %.1f, %.1f), Semi-dims = (%.1f, %.1f, %.1f)\n", ...
        X_0, Y_0, Z_0, a, b, c);

fprintf("Interior Laplacian = -4πGρ = %.6f\n", -4*pi*G*rho);
fprintf("--------------------------------------------------------------------------\n");

% Print table header
fprintf("%-28s %-40s %-12s %-12s %-12s %-12s\n", ...
    "Test Point (Abs)", "Description", "Vxx", "Vyy", "Vzz", "Laplacian");

% Print rows
for k = 1:numTests
    fprintf("%-28s %-40s %12.8f %12.8f %12.8f %12.8f\n", ...
        results(k).AbsPoint, ...
        results(k).Description, ...
        results(k).Vxx, results(k).Vyy, results(k).Vzz, results(k).Laplacian);
end
