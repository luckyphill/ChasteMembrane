function x=tear_off(l,n,del,s,cl,x0,y0)
	% A function returning the resting length X that is the x coordinate distance from the membrane wall of a single isolated epithelial cell
	% l - natural spring length
	% n - number of membrane cells upwards and downwards  - total number will be 2n+1
	% del - the spacing between each membrane cell
	% x0 - a guess at the resting position - this is used as a proxy for starting position, it seems to work correctly

	% Need to make a 'switch-finding' method to find the point where a zero no longer exists
	pfx_upper = 100;
	pfx_lower = 0.01;
	% The aim is to minise the difference between upper and lower, to a specified tolerance
	difference = 10;
	tol = 1e-6;
	% Also provide an iteration limit in case of non-convergence
	iterations = 0;
	it_limit = 100;


	options = optimset('Display','off');
	[~,~,exit_flag_upper] = fsolve(@(X)net_force(X,l,n,del,s,cl,pfx_upper,0), [x0,y0], options);
	[~,~,exit_flag_lower] = fsolve(@(X)net_force(X,l,n,del,s,cl,pfx_lower,0), [x0,y0], options);

	% Make sure the upper bound is large enough
	while exit_flag_upper == 1
		pfx_upper = 2 * pfx_upper;
		[~,~,exit_flag_upper] = fsolve(@(X)net_force(X,l,n,del,s,cl,pfx_upper,0), [x0,y0], options);
	end

	% Make sure the lower bound is large enough
	while exit_flag_lower ~= 1
		pfx_lower = pfx_lower / 2;
		[~,~,exit_flag_lower] = fsolve(@(X)net_force(X,l,n,del,s,cl,pfx_lower,0), [x0,y0], options);
	end


	
	while difference > tol && iterations < it_limit
		pfx_new = mean([pfx_upper,pfx_lower]);
		[~,~,exit_flag_new] = fsolve(@(X)net_force(X,l,n,del,s,cl,pfx_new,0), [x0,y0], options);
		if exit_flag_new ~=1
			pfx_upper = pfx_new;
		end
		if exit_flag_new ==1
			pfx_lower = pfx_new;
		end

		if pfx_new ~= pfx_lower && pfx_new ~= pfx_upper
			fprintf('Something has gone wrong in the iteration')
			break;
		end

		difference = abs(pfx_upper-pfx_lower);
		iterations = iterations + 1;
	end

	x = pfx_lower;

end