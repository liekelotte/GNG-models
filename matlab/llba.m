%%%%%%%%%%%%%%%%%
%%% Written by Quentin J. M. Huys, UCL, London 2011
%%% Reference:
%%% Guitart-Masip M, Quentin JM, Fuentemilla LL, Dayan P, Duzel E, Dolan RJ (2012)
%%% Go and no-go learning in reward and punishment: Interaction between affect and effect NeuroImage doi:10.1016/j.neuroimage.2012.04.024

function [l] = llba(x,a,r,s,Z,doprior)

beta = exp(x(1));            % sensitivity to reward          
alpha = 1./(1+exp(-x(2)));     % learning rate

if doprior
	l = -1/2 * (x-Z.mu)'*Z.nui*(x-Z.mu) -  1/2*log(1/det(Z.nui/(2*pi))); %
else
	l=0;
end


V=zeros(1,4); 
Q=zeros(2,4); 


for t=1:length(a)
	er = beta * r(t);

	q = Q(:,s(t)); 

	l0 = max(q);
	la = q(a(t)) - l0 - log(sum(exp(q-l0)));
	l = l + la;


	Q(a(t),s(t)) = Q(a(t),s(t)) + alpha * (er - Q(a(t),s(t)));  

end
l  = -l ;

