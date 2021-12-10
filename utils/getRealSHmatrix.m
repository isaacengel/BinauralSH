function Y = getRealSHmatrix(N,az,el)
% Get matrix with real spherical harmonics.
% INPUT:
%   N = SH order
%   az = azimuth (ndirs x 1) in rad
%   el = elevation (ndirs x 1) in rad (0=top, pi/2=front)
%   
% OUTPUT:
%   Y = real SH matrix ((N+1)^2 x ndirs)
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% August 2021

az = az(:).'; % force row
el = el(:).';

nsh = (N+1)^2;
ndirs = size(az,2);

Y = zeros(nsh,ndirs);

for n=0:N
    for m=-n:n
        ind = n^2 + n + m + 1; % ACN ordering
        % Elevation term
        Pn = legendre(n,cos(el.'));
        % Azimuth term
        if m<0
            ym = sqrt(2)*sin(abs(m)*az);
        elseif m>0
            ym = sqrt(2)*cos(m*az);
        else
            ym = ones(1,ndirs);
        end
        % Normalisation term (N3D)
        Norm = sqrt( (2*n+1)/(4*pi) * factorial(n-abs(m))/factorial(n+abs(m)));
        % Final formula
        Y(ind,:) = (-1)^m * Norm .* Pn(abs(m)+1,:) .* ym;
    end
end

