function [Xc] = OptimalShrinkage(X, sigma)

    N = size(X, 1) ; 
    p = size(X,2) ;
	beta = p ./ N ;
    
	[u,l,v] = svd(X./sqrt(N)./sigma); 
    J = zeros(size(l)); 
	y = diag(l) ; eta = sqrt( (y.^2-beta-1).^2 - 4*beta) ./ y ;
	tmp = find(y<=1+sqrt(beta)) ; eta(tmp) = 0 ;
	tmp = min(size(J)) ; J(1:tmp,1:tmp) = diag(eta) ; 
	Xc = (sigma*sqrt(N))*u*J*v' ; 
    
end