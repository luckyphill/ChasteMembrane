
steps = 100;
xf = linspace(0,5,steps);

%F(steps) = struct('cdata',[],'colormap',[]);
for j = 1:steps
	h = grad_map(1,10,.5,10,2,0,xf(j));
    F(j) = getframe(h);
end

fig = figure;
movie(fig,F,2)