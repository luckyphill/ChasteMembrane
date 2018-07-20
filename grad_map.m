function varargout = grad_map(l,n,del,s,cl,pfx,pfy)

	res = 1000;
	x = linspace(.1,cl-.1,res);
	y = linspace(-1,1,res);

	grad_plot = zeros(100);

	for i =1:res
		for j= 1:res
			f = net_force([x(i), y(j)],l,n,del,s,cl,pfx,pfy);
			grad_plot(i,j) = (f(1)^2 + f(2)^2)^(1/2);
		end
	end
	h = figure();
	imagesc(rot90(grad_plot),[0 30]);
	colorbar;
	varargout{1} = h;

end