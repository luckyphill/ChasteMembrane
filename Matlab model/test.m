clear all
close all

% A script that solves the 1D layer of proliferating cells model
% The positions of the cells are stored in the matrix p.x
% At each time step a new row (column?) of the matrix is added
% When a new cell is added, the matrix shifts its columns (rows?) across to
% accommodate the new cell.
% Cells are killed above a certain limit - p.top
% When a cell is killed is no longer passed into the force calculation but
% it still takes space in the matrix.
% In order to plot things correctly, the positions of new cells before they
% divide and the positions of old cells after they are killed are stored as
% nans.
% 

p.n = 20; % the initial number of cells
p.n_dead = 0; % number of cells that have died

p.t_end = 200;
p.dt = 0.02;

p.x = 0:p.n-1; % intial positions
p.v = zeros(size(p.x)); % initial velocities for plotting only

p.ages = 10 * rand(1,p.n); % randomly assign ages at the start
p.divide_age = get_a_divide_age(p.n); % randomly assign an age when division occurs
p.divide_age(1) = p.t_end + 14; % a quick hack to stop bottom cell dividing

p.division_spring_length = 0.1; % after a cell divides, the new cells will be this far apart
p.growth_time = 1.0; % time it takes for newly divided cells to grow to normal disatance apart
p.cut_out_height = 15; % the height where proliferation stops

p.ci = true;
p.ci_fraction = 0.88; % the compression on the cell that induces contact inhibition as a fraction of the free volume

p.l = 1; % The natural spring length of the connection between two mature cells
p.k = 20; % The spring constant
p.damping = 1.0; % The damping constant
p.top = 20; % The position of the top of the wall

assert(p.top>=p.cut_out_height);

p.t_start = 0; % starting time

p.t = p.t_start;

while p.t < p.t_end
    
    % Next time step, and age the cells
    p.t = p.t + p.dt;
    p.ages = p.ages + p.dt;
       
    % Refresh the indices of alive and dead cells
    alive = 1:p.n;
    dead = p.n+1:p.n+p.n_dead;
    
    % Force calculations only for the live cells
    f = force(p.x(end,alive),p);
    
    % New positions for live cells
    p.x(end+1,alive) = p.x(end,alive) + f * p.dt/p.damping;
    p.x(end,dead) = nan(1,p.n_dead); % bookkeeping for the dead cells
    
    p.v(end+1,alive) =  f/p.damping;
    p.v(end,dead) = nan(1,p.n_dead); % bookkeeping for the dead cells (not end+1 because previous command created the new row)
    
       
    cells_to_divide = []; % reset the list of cells to divide
    temp = 1:p.n; % used to get the indices
    proliferative_zone = temp(p.x(end,:)<p.cut_out_height); %determines cells that can proliferate
    cells_to_divide = temp(p.ages(proliferative_zone) > p.divide_age(proliferative_zone)); % determines cells ready to divide

    % Process the cells ready to divide
    if length(cells_to_divide) > 0
        p = divide_cells(cells_to_divide,p);
    end
    
    % Kill any cells past the top of the crypt
    p = sloughing(p);
    
    
    
end

plot_cells(p)
