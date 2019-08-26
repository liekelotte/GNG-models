%%%%%%%%%%%%%%%%%
%%% Written by Quentin J. M. Huys, UCL, London 2011
%%% Reference:
%%% Guitart-Masip M, Quentin JM, Fuentemilla LL, Dayan P, Duzel E, Dolan RJ (2012)
%%% Go and no-go learning in reward and punishment: Interaction between affect and effect NeuroImage doi:10.1016/j.neuroimage.2012.04.024

function [l] = ll2baxb(x, a, r, s, Z,doprior)
beta 	  = exp(x(1:2));            % sensitivity to reward          
alfa 	  = 1./(1+exp(-x(3)));     % learning rate             % 'pavlovian' parameter. Weigth of Vcue into Qgo
g       = 1/(1+exp(-x(4)));
bias	  = x(5);

if doprior
	 lp = -1/2 * (x-Z.mu)'*Z.nui*(x-Z.mu) -  1/2*log(1/det(Z.nui/(2*pi))); %
else
	lp = 0;
end

% initialize 
l=0;

Q=zeros(2,4); 

for t=1:length(a)
	rho = sum(s(t)==[1 3]);

	er = beta(2-rho) * r(t);

	q = Q(:,s(t)); 
	q(1) = q(1) + bias;   

	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;
	l = l + log(pg(a(t)));

	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa * (er - Q(a(t),s(t)));  
end
l  = -l  - sum(lp);



