function [F_R, theta,F_Ry,F_Rx] = forcesresultant(m1,m2)


% see http://secs.oakland.edu/~latcha/EGR280/StatDynMatlab.pdf

% limit impact of very small (noisy) coupling
if abs(m1) < .001
    m1 = 0;
end
if abs(m2) < .001
    m2 = 0;
end

F_x1 = m2;
F_y1 = m2;
F_x2 = -m1;
F_y2 = m1;


% Sum the y-components of the two forces to determine the y-component of the resultant force
F_Ry = F_y1 + F_y2;
% Sum the x-components of the two forces to determine the x-component of the resultant force
F_Rx = F_x1 + F_x2;

% Calculate the magnitude of the resultant force
F_R  = sqrt( F_Rx^2+F_Ry^2);

% Calculate the angle of the resultant force% (in degrees from -90)
u = [0 -1];% -90 deg
v = [F_Rx F_Ry];
x1 = u(1);
x2 = v(1);
y1 = u(2);
y2 = v(2);
CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
theta = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);

