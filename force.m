function f = force(x,p)

% takes a vector of cell centre positions
% assuming a 1D column of cells
% returns the net force on each cell by using the positions to calculate
% the spring forces


j=2:p.n;
lengths = x(j) - x(j-1); % the lengths of each spring

age_comparison = p.ages(j) - p.ages(j-1);
Dx = lengths - p.l; % vector of displacements from natural spring length

%Need to adjust natural spring length for newly divided cells
for i=1:p.n-1
    if (age_comparison(i) == 0.0) && (p.ages(i+1) < p.growth_time) % This has a bug if two adjacent cells divide at the same time
        Dx(i) = lengths(i) - (p.division_spring_length + (p.l - p.division_spring_length) * p.ages(i+1) / p.growth_time);
    end
end
    


f_s = spring_force(Dx,p.k,p.l); % -ve force means spring is compressed, +ve force means spring is in tension

% tension pulls below cell up and above cell down
% compression pushes below cell down and above cell up
% tension spring stored as +ve force
% compression spring stored as -ve force

m = 2:p.n-1;

f(1) = 0; % first cell is fixed in space
f(m) = f_s(m) - f_s(m-1);
f(p.n) = -f_s(p.n-1);


end

function f = spring_force(dx,k,l)
    % function to return spring force
    % this way allows for easy modification of force rule
    % a negative force means the spring is in tension
    
    %linear spring
    %f = k*dx;
    
    for i=1:length(dx)
        if dx(i)>0
            f(i) = k * dx(i) * exp(-1.8 * dx(i)/l);
        else
            f(i) = k * l * log(1.0 + dx(i)/l);
        end
    end
end