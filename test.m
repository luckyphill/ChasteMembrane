clear all
close all
x = [0,1,2,3,4,5,6,7,8,9,10];
v = zeros(size(x));
dt = 0.1;
l=0.9;
k=5;
t_end = 100;
t = 0;

while t < t_end
    f = force(x(end,:),l,k);
    [dx, v_new] = displacement(f,v(end,:),dt);
    x(end+1,:) = x(end,:) + dx;
    v(end+1,:) = v_new;
    t = t + dt;
end