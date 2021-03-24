function plotSph(az,el,X)
% Plot data in vector X defined as a function of azimuth and elevation 
% (e.g. to draw a loudness map of an HRTF), interpolating values as needed.
%
% INPUT:
%   az = vector of azimuths in rad (0=front, pi/2=left)
%   el = vector of elevations in rad (0=top, pi/2=front)
%   X = vector of data, same length as az and el
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021
  
x = az(:)*180/pi ; y = 90-el(:)*180/pi ; z = X(:);
% concatenate a second "lap" with negative azimuths to help the plot
x = [x;x-360]; 
y = [y;y];
z = [z;z];
% triangulate and plot as a surface
dt = delaunayTriangulation(x,y) ;
tri = dt.ConnectivityList ;
xi = dt.Points(:,1) ; 
yi = dt.Points(:,2) ; 
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi) ;
trisurf(tri,xi,yi,zi) 
set(gca,'XDir','reverse','XLim',[-180 180],'YLim',[-88 88])
view(2)
shading interp