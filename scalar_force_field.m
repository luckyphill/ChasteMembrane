
steps = 100;
xf = linspace(0,10,steps);

%F(steps) = struct('cdata',[],'colormap',[]);
for j = 1:steps
	h = grad_map(1,50,.1,2,2,xf(j),0);
    F(j) = getframe(h);
end

fig = figure;
movie(fig,F,2)