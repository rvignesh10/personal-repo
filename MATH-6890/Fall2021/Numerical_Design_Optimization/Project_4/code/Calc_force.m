%% Calculate force distribution along length of spar

function force = Calc_force(x,Mass,L)
% Input - x (length discretization), Mass, Length of spar;
% Output - force at each node

% force is linear in nature - cX 
% integral(cx dx)=2.5*Weight/2 - c = 2.5*Weight/L^2

g = 9.81; % m/s^2 
c = 2.5*Mass*g/L^2;
force = c*x;
force = flip(force);
end