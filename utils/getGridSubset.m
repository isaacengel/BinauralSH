function [gridout,logic_ind,dist_vec] = getGridSubset(gridin,gridreq,doplot)
% Return a subset of points within a grid that are closest to a set of
% requested points.
%
% INPUT:
%   gridin = original grid (npoints x 2). First column is azimuth in rad,
%       second column is elevation in rad (0=top)
%   gridreq = matrix of requested points, same format as gridin
%   doplot = do some plots
%
% OUTPUT:
%   gridout = output grid, same format as gridin
%   logic_ind = indices of the output grid within the original one
%   dist_vec = vector with distances from each output point to the
%       corresponding requested one 
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

azin = gridin(:,1);
elin = pi/2-gridin(:,2); 
azreq = gridreq(:,1);
elreq = pi/2-gridreq(:,2); 

N = size(gridreq,1);
azout = zeros(N,1);
elout = zeros(N,1);
logic_ind = false(size(gridin,1),1);
gridout = zeros(N,2);
dist_vec = zeros(N,1);

for i=1:N
    dist = acos(sin(elin).*sin(elreq(i))+cos(elin).*cos(elreq(i)).*cos(azin-azreq(i)));
    [dist_vec(i),ind] = min(dist);
    azout(i) = azin(ind);
    elout(i) = elin(ind);
    gridout(i,:) = [azout(i),pi/2-elout(i)];
    logic_ind(ind) = 1;
end

if doplot
    figure
    [Xm1,Ym1,Zm1]=sph2cart(azreq,elreq,1.01);
    [Xm2,Ym2,Zm2]=sph2cart(azout,elout,1.02);
    colormap Gray;
    plot3(Xm1,Ym1,Zm1,'marker','o','markerfacecolor','r','color','r','linestyle','none','MarkerSize',6)
    hold on
    plot3(Xm2,Ym2,Zm2,'marker','d','markerfacecolor','g','color','g','linestyle','none','MarkerSize',6)
    axis off; hold on; grid off;
    sphere;
    axis equal; rotate3d on;
    light;
    alpha(.8);
    lighting phong;
%     camzoom(1.4);
    hold off;
    title('Requested (red) and returned (green) grids')
end
