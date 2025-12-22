%% ============================================================
% Example 03 — Gravitational Acceleration Field (z = 0 Plane)
% ============================================================

clear; clc;

%% Add utils
addpath('/Volumes/Dunendran/Programs/MATLAB/Cuboid_Extension/utils_MATLAB');

%% ---------------------------------------------------------
% CUBE PARAMETERS
% ---------------------------------------------------------
G   = 1;
rho = 1;

X0 = 0;
Y0 = 0;
Z0 = 0;

L = 2.0;
B = 2.0;
D = 2.0;

% Semi-dimensions
a = L/2;
b = B/2;
c = D/2;

%% ---------------------------------------------------------
% GRID (z = 0 PLANE)
% ---------------------------------------------------------
x_vals = linspace(-2, 2, 13);
y_vals = linspace(-2, 2, 13);
[X, Y] = meshgrid(x_vals, y_vals);

U = zeros(size(X));
V = zeros(size(Y));

%% ---------------------------------------------------------
% ANALYTICAL ACCELERATION FIELD
% ---------------------------------------------------------
for i = 1:size(X,1)
    for j = 1:size(X,2)

        gx = GravCuboidAcceleration.accel_gx(X(i,j), Y(i,j), 0, a,b,c, rho, G);
        gy = GravCuboidAcceleration.accel_gy(X(i,j), Y(i,j), 0, a,b,c, rho, G);

        U(i,j) = gx;
        V(i,j) = gy;
    end
end

%% ---------------------------------------------------------
% scale=65
% ---------------------------------------------------------
U_plot = U / 65;
V_plot = V / 65;

%% ---------------------------------------------------------
% FIGURE — figsize=(8,7)
% ---------------------------------------------------------
fig = figure('Color','w','Units','inches','Position',[1 1 8 7]);
ax = axes('Parent',fig); hold(ax,'on');

%% ---------------------------------------------------------
% QUIVER 
% ---------------------------------------------------------
q = quiver(ax, X, Y, U_plot, V_plot, ...
    'Color','red', ...
    'LineWidth', 0.003 * 72);   % width=0.003 → convert to points

% Remove quiver from legend (fixes unwanted "data")
q.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Improve arrowhead to emulate pivot='middle'
q.MaxHeadSize = 0.5;

%% ---------------------------------------------------------
% CUBOID BOUNDARY
% ---------------------------------------------------------
x_min = X0 - L/2;
x_max = X0 + L/2;
y_min = Y0 - B/2;
y_max = Y0 + B/2;

plot(ax, ...
    [x_min x_max x_max x_min x_min], ...
    [y_min y_min y_max y_max y_min], ...
    'k-', 'LineWidth', 1, 'DisplayName','Cube Boundary');

%% ---------------------------------------------------------
% LEGEND 
% ---------------------------------------------------------
lg = legend(ax, 'Location','southwest', 'Box','off');
lg.TextColor = 'black';
lg.FontSize  = 16;

%% ---------------------------------------------------------
% AXES STYLE
% ---------------------------------------------------------
title(ax, "Gravitational Acceleration in the Plane z = 0", ...
      'FontSize', 17, 'Color','black');

xlabel(ax, 'x', 'FontSize',16,'Color','black');
ylabel(ax, 'y', 'FontSize',16,'Color','black');

ax.XLim = [-2.0 2.0];
ax.YLim = [-2.0 2.0];
ax.XTick = -2:1:2;
ax.YTick = -2:1:2;

ax.FontSize = 16;
ax.XColor = 'black';
ax.YColor = 'black';

axis(ax, 'equal');
grid(ax, 'off');
box(ax, 'on');

%% ---------------------------------------------------------
% REMOVE WHITESPACE
% ---------------------------------------------------------
outer = ax.OuterPosition;
ti = ax.TightInset;
ax.Position = [
    outer(1) + ti(1), ...
    outer(2) + ti(2), ...
    outer(3) - ti(1) - ti(3), ...
    outer(4) - ti(2) - ti(4)];

%% ---------------------------------------------------------
% SAVE
% ---------------------------------------------------------
exportgraphics(fig, "Gravitational_Acceleration_quiver_z0_matlab.png", ...
    'Resolution',300);

disp("Saved: Gravitational_Acceleration_quiver_z0_matlab.png");
