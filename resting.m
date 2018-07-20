function x=resting(l,n,del,s,cl, x0,pfx,pfy)
	% A function returning the resting length X that is the x coordinate distance from the membrane wall of a single isolated epithelial cell
	% l - natural spring length
	% n - number of membrane cells upwards and downwards  - total number will be 2n+1
	% del - the spacing between each membrane cell
	% x0 - a guess at the resting position - this is used as a proxy for starting position, it seems to work correctly

	options = optimset('Display','off');
	[a,~,exit_flag] = fsolve(@(X)net_force(X,l,n,del,s,cl,pfx,pfy), x0, options);
	x = a;


end