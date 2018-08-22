function p = divide_cells(dividing_cells,p)
    
    % takes in a vector of cells that are ready to divide
    % inserts the new cell into the matrices/vectors next to the parent
    % cell
    
    % newly divided cell appears ABOVE parent cell
    
    for i = 1:length(dividing_cells)
        index = dividing_cells(i);
        
        % deal with age tracking
        p.ages = [p.ages(1:index-1), 0,0, p.ages(index+1:end)];
        p.divide_age = [p.divide_age(1:index-1), get_a_divide_age(2), p.divide_age(index+1:end)];
        
        % insert a new row into the vector of positions and velocities
        p.x  = [p.x(:,1:index), nan(size(p.x(:,1))), p.x(:,index+1:end)];
        p.v  = [p.v(:,1:index), nan(size(p.x(:,1))), p.v(:,index+1:end)];
        
        predivision = p.x(end,index);
        p.x(end,index) = predivision - 0.5 * p.division_spring_length;
        p.x(end,index+1) = predivision + 0.5 * p.division_spring_length;
        
        % This is a bit of a hack, but should work
        % Once the new cell has been inserted, the divide indices will be
        % wrong
        % They should be in order of the indices lowest to highest. If this
        % is the case then the newly inserted cell will always be below all the
        % cells that are yet to divide, so we can adjust the indices by
        % adding 1 to them all
        
        p.n = p.n + 1; % increase number of cells
        
        dividing_cells = dividing_cells + 1; % correct indexing error
        
        
    end


end