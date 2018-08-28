function p = divide_cells(dividing_cells,p)
    
    % takes in a vector of cells that are ready to divide
    % inserts the new cell into the matrices/vectors next to the parent
    % cell
    
    % newly divided cell appears ABOVE parent cell
    
    pause_time = 4;
    
    if p.ci
        vol = get_cell_volume(p);
        comp_cells = [];
        for i = 1:length(dividing_cells)
            if (vol(i) < p.ci_fraction * 2 * p.l)
                comp_cell = dividing_cells(i);
                % remove this cell from consideration
                % there are three ways to do this
                
                
                % 1. skip division this cycle
                %p.divide_age(comp_cell) = p.divide_age(comp_cell) + get_a_divide_age();
                
                % 2. pause the cell for x hours then trying again
                p.divide_age(comp_cell) = p.divide_age(comp_cell) + pause_time;
                
                % 3. pause the cell until the compression is gone
                %p.divide_age(comp_cell) = p.divide_age(comp_cell) + p.dt;
                
                comp_cells = [comp_cells, i];
                
            end
        end
        dividing_cells(comp_cells) = [];
    end
    
    
    for i = 1:length(dividing_cells)
        
        cell_d = dividing_cells(i); % Index of the dividing cell
        cell_n = cell_d + 1; % Index of the new cell
        above_n = cell_n:(p.n+p.n_dead); % The cells above the new cell
        below_d = 1:(cell_d-1); % The cells below the dividing cell
        below_n = 1:cell_d; % The cells below the new cell (includes dividing cell)
        
        fprintf('%d\n',length(p.x(:,cell_d)));
        % deal with age tracking
        p.ages = [p.ages(below_d), 0,0, p.ages(above_n)];
        p.divide_age = [p.divide_age(below_d), get_a_divide_age(2), p.divide_age(above_n)];
        
        % insert a new row into the vector of positions and velocities
        p.x  = [p.x(:,below_n), nan(size(p.x(:,1))), p.x(:,above_n)];
        p.v  = [p.v(:,below_n), nan(size(p.x(:,1))), p.v(:,above_n)];
        
        predivision = p.x(end,cell_d);
        fprintf('%d\n',length(p.x(:,cell_d)));
        p.x(end,cell_d) = predivision - 0.5 * p.division_spring_length;
        p.x(end,cell_n) = predivision + 0.5 * p.division_spring_length;
        
        % This is a bit of a hack, but should work
        % Once the new cell has been inserted, the divide indices will be
        % wrong
        % They should be in order of the indices lowest to highest. If this
        % is the case then the newly inserted cell will always be below all the
        % cells that are yet to divide, so we can adjust the indices by
        % adding 1 to them all
        
        p.n = p.n + 1; % increase number of cells
        
        dividing_cells = dividing_cells + 1; % correct indexing error
        fprintf('Cell %d divided at t = %.2f\n', cell_d,p.t);
        
    end


end