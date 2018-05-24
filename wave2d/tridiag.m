function u = tridiag(a,b,c,r,N)
% Function tridiag: 
%    Inverts a tridiagonal system whose lower, main and upper diagonals
%    are respectively given by the vectors a, b and c. r is the right-hand
%    side, and N is the size of the system. The result is placed in u.
%    (Adapted from Numerical Recipes, Press et al. 1992)
if (b(1)==0) 
    fprintf(1,'Reorder the equations for the tridiagonal solver...') 
    pause
end
beta = b(1) ;
u(1) = r(1)/beta ;
% Start the decomposition and forward substitution
for j = 2:N
    gamma(j) = c(j-1)/beta ;
    beta = b(j)-a(j)*gamma(j) ;
    if (beta==0)
        fprintf(1,'The tridiagonal solver failed...') 
        pause
    end
    u(j) = (r(j)-a(j)*u(j-1))/beta ;
end
% Perform the backsubstitution
for j = 1:(N-1)
    k = N-j ;
    u(k) = u(k) - gamma(k+1)*u(k+1) ;
end