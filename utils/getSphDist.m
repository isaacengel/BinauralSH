function dist = getSphDist(azA,elA,azB,elB)
% Great circle distance between two points in spherical coordinates (rad)
elA = pi/2 - elA; % elevation is expected as (pi/2=front) but the formula below assumes (0=front)
elB = pi/2 - elB;
dist = acos(sin(elA).*sin(elB)+cos(elA).*cos(elB).*cos(azA-azB));


