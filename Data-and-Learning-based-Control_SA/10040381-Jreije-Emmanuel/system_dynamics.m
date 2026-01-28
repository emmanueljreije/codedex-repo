function x_next = system_dynamics(x, u)
A = [1, -18.4542,-3.3219, 2.2083, 0.1722;
     0, 0.6807, 0.2618, 0, 0;
     0, 0, 0.6807, 0, 0;
     0, 0, 0, 0.8519, 0.1365;
     0, 0, 0, 0, 0.8519];
B = [-0.4539, 0.0094;
    0.0575, 0;
    0.3193, 0;
    0, 0.01155;
    0, 0.1481];

if size(x,1) ~= 5
    fprintf('The state input does not have the correct dimension! Stopping program \n');
    return
elseif size(x,2) ~=1
    fprintf('The state input does not have the correct dimension! Stopping program \n');
    return
end
 
if size(u,1) ~= 2
    fprintf('The control input does not have the correct dimension! Stopping program \n');
    return
elseif size(u,2) ~= 1
    fprintf('The control input does not have the correct dimension! Stopping program \n');
    return
end

% w = 0.1*rand;
% w = 0;
x_next = A*x + B*u;% + w;
end

