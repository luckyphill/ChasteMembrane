function plot_cells(x,r, sampling_multiple)
% takes a n by t_steps matrix which contains the positions of n cells at
% each of t_steps time steps
% for each cell it draws a circle of specified radius r
% it then draws the circles from the matrix at the time steps
% dt*sampling_multiple

[t_steps, n] = size(x);

figure
cla
X = zeros(n,1);
radii = (r/2) * ones(n,1);

xlim([-n*r/2 n*r/2])
ylim([-1 (n+1)*r])
axis square
for t=1:sampling_multiple:t_steps
    cla
    centres = [X x(t,:)'];
    viscircles(centres,radii);
    pause(0.02);

end