function plot_cells(p)
% takes a struct containing all the data about the cell population

[t_steps, n] = size(p.x);

figure('pos',[1 165 960 892])
cla
X = zeros(n,1);
radii = (p.l/2) * ones(n,1);

lims = ((p.top+1)*p.l +1)/2;

xlim([-lims lims])
ylim([-1 (p.top+1)*p.l])
axis square
for t=1:t_steps
    cla
    centres = [X p.x(t,:)'];
    viscircles(centres,radii);
    title((t-1)*p.dt)
    pause(p.dt);

end