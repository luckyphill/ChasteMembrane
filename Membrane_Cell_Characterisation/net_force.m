function f_net = net_force(X,l,n,del,s,cl,push_force_x,push_force_y)
	% A function returning the net force on an isolated epithelial cell
	% X is a vector containing the position of the epithelial cell X = (x,y)

	y = zeros(2*n+1,1);
	for i=1:n
		y(i+1) = +i*del;
		y(n + i +1) =  - i*del;
	end

	d = ((y-X(2)).^2 + X(1)^2).^(1/2);

	f = zeros(2*n+1,1);
	for i = 1:2*n+1
		f(i) = force(d(i), l, s, cl);
	end
	f_net = zeros(2,1);
	f_net(1) = -sum(X(1).*f./d) + push_force_x; % needs negative in front of sum since attraction is considered +ve and that will happen in the -ve x direction
	f_net(2) = sum(-(y-X(2)).*f./d) + push_force_y;

end


function f = force(d,l,s, cl)
	% A function calculating the signed force given a natural length (l) and a current length (d)
	% s is the spring stiffness
	% b is a scaling parameter determining the attraction curve
	% cl is the maximum sensing distance, so beyond this cells do not interact. cl = cutoff length

	b = 1.8;
	f = 0;
	if d <= cl
		delta = d - l;
		if delta > 0
			f = s * delta * exp(-b * delta/l);
		else
			f = s * l * log(1 + delta/l);
		end
	end

end