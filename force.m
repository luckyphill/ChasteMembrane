function f = force(x,l,k)

% takes a vector of cell centre positions
% assuming a 1D column of cells
% returns the net force on each cell by using the positions to calculate
% the spring forces

n = length(x);
% n cells means n-1 springs

j=2:n;
lengths = x(j) - x(j-1); % the lengths of each spring

Dx = lengths - l; % vector of displacements from natural spring length

f_s = spring_force(Dx,k); % -ve force means spring is compressed, +ve force means spring is in tension

% tension pulls below cell up and above cell down
% compression pushes below cell down and above cell up
% tension spring stored as +ve force
% compression spring stored as -ve force

m = 2:n-1;

f(1) = 0; % first cell is fixed in space
f(m) = f_s(m) - f_s(m-1);
f(n) = -f_s(n-1);



end

function f = spring_force(dx,k)
    % function to return spring force
    % this way allows for easy modification of force rule
    % a negative force means the spring is in tension
    
    %linear spring
    f = k*dx;
    
end