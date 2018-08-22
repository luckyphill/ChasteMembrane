clear all
close all

p.n = 11; % store the number of cells
p.n_dead = 0; % number of cells that have died

p.t_end = 20;
p.dt = 0.01;

p.x = 0:p.n-1; % intial positions
p.next_index = p.n+1; % contains the new index of a newly born cell
p.nodes = 1:p.n; % vector containing the node indices in order they appear in x

p.ages = 8 * rand(1,p.n); % randomly assign ages at the start
p.divide_age = get_a_divide_age(p.n); % randomly assign an age when division occurs
p.divide_age(1) = p.t_end; % a quick hack to stop issues

p.division_spring_length = 0.05; % after a cell divides, the new cells will be this far apart
p.growth_time = 1.0; % time it takes for newly divided cells to grow to normal disatance apart

p.l = 1;
p.k = 15;
p.damping = 1.0;
p.top = 15;

t = 0;

while t < p.t_end
    % function to kill cells past the top
    p = sloughing(p);
    alive = 1:p.n;
    dead = p.n+1:p.n+p.n_dead;
    
    f = force(p.x(end,alive),p);
    
    p.x(end+1,alive) = p.x(end,alive) + f * p.dt/p.damping;
    p.x(end,dead) = nan(1,p.n_dead); % NaNs keep the matrix full, but stop plotting
    
    t = t + p.dt;
    p.ages = p.ages + p.dt;
    
    cells_to_divide = []; % reset the list of cells to divide
    temp = 1:p.n;
    cells_to_divide = temp(p.ages(alive) > p.divide_age(alive));
    
    p = divide_cells(cells_to_divide,p);
    % process the cells that are ready to divide
    
end
figure('pos',[1 165 960 892])
plot(0:p.dt:p.t_end,p.x)
plot_cells(p)
