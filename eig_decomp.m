function [V,D] = eig_decomp(X, method, niter)

%function that does eigenvalue decomposition with the QR algorithm
%V - eigenvectors, D = diagonal matrix with eigenvalues
%X - input matrix (must be square)
%niter = number of iterations
%method - algorithm to be used :
%gram_schmidt - gram schmidt QR O(n^3)
%hessenberg - hessenberg QR O(n^2)
%tridiagonal - symmetric tridiagonal QR with implicit Wilkinson shift

%notes - in general, U returns the Schur vectors, not the eigenvectors
%but for normal matrices (that includes hermitian matrices), the schur
%vectors are the eigenvectors

if(nargin == 1)
    method = 'gram_schmidt';
elseif(nargin == 2)
    niter = 100;
end

A = X;
n = length(A);

if (strcmp(method,'gram_schmidt')) 
    U = eye(n);
    for i = 1:niter
        [Q,R] = gram_schmidt(A);
        A = R*Q;
        U = U*Q;
    end
    D = diag(diag(A));
    V = U;
    
elseif(strcmp(method,'hess'))
        [H,U] = hessenberg(A);
        c = zeros(n-1,1);
        s = zeros(n-1,1);
        for i = 1:niter
            for k = 1:n-1
                %H = G(k,k+1,v)*H
                [c(k),s(k)] = givens(H(k,k),H(k+1,k));
                rot_mat = [c(k), -s(k);s(k), c(k)];
                H(k:k+1,k:n) = rot_mat * H(k:k+1,k:n);
            end
            for k = 1:n-1
                %apply givens rotation from right
                rot_mat = [c(k), s(k);-s(k), c(k)];
                H(1:k+1,k:k+1) = H(1:k+1,k:k+1)*rot_mat;
                U(:,k:k+1) = U(:,k:k+1)*rot_mat;
            end
        end
        V = U;
        D = diag(diag(H));
        
elseif (strcmp(method,'tridiag'))
    
    [T,U] = hessenberg(A,'sym');
    %main diagonal elements of A
%    a = diag(T);
%    %off-diagonal elements
%    b = [0;diag(T(2:n,1:n-1))];
%    iter = 0;
%    eps = 10^-6;
    
%     %m reduced in convergence check
%     m = n;
%     while(m > 1)
%         iter = iter+1;
%         %compute Wilkinson shift
%         d = (a(m-1) - a(m))/2;
%         if d == 0
%             s = a(m)-abs(b(m));
%         else
%             s = a(m) - (b(m)*b(m))/(d + sign(d)*safeDistance(d,b(m)));
%         end
%         %implicit QR step begins
%         x = a(1)-s;
%         y = b(2);
%         for k = 1:m-1
%             if m > 2
%                 dist = safeDistance(x,y);
%                 c = x/dist;
%                 s = -y/dist;
%                 %[c,s] = givens(x,y);
%             else
%                 alpha = (a(1) - a(2))/b(2);
%                 denom = safeDistance(1,alpha);
%                 c = alpha/denom;
%                 s = -1/denom;     
%             end
%             w = c*x - s*y;
%             d = a(k) - a(k+1);
%             z = (2*c*b(k+1) + d*s)*s;
%             a(k) = a(k) - z;
%             a(k+1) = a(k+1) + z;
%             b(k+1) = d*c*s + (c*c - s*s)*b(k+1);
%             x = b(k+1);
%             if k>1
%                 b(k) = w;
%             elseif k < m-1
%                 y = -s*b(k+2);
%                 b(k+2) = c*b(k+2);
%             end
%             U(:,k:k+1) = U(:,k:k+1)*[c,s;-s,c];
%         %end implicit QR step
%         end 
%         
%         %check for convergence
%         if abs(b(m)) < eps*(abs(a(m-1)) + abs(a(m))) || (iter >= niter)
%             if(iter >= niter)
%                 warning('ImplicitTriQR:instability',['ImplicitTriQR() did not converge for m = ' num2str(m)]);
%             end
%             m = m-1;
%             iter = 0;
%         end 
%     end
    
% Perform implicit QR algorithm
d = diag(T);
od = [0;diag(T(2:n,1:n-1))];
m = n;
iter = 0;
eps = 10^-6;

%this part of the code taken from Linear Algebra Package by Brian Moore
while (m > 1)
    iter = iter + 1;
    g = (d(m-1) - d(m)) / 2;
    if (g == 0)
        s = d(m) - abs(od(m));
    else
        s = d(m) - od(m) * od(m) / (g + sign(g) * SafeDistance(g,od(m)));
    end
    x = d(1) - s;
    y = od(2);
    for k = 1:(m-1)
        if (m > 2)
            xydist = SafeDistance(x,y);
            c = x / xydist;
            s = -y / xydist;
        else
            alpha = (d(1) - d(2))/od(2);
            denom = SafeDistance(1,alpha);
            c = alpha / denom;
            s = -1 / denom;
        end
        w = c * x - s * y;
        g = d(k) - d(k+1);
        z = (2 * c * od(k+1) + g * s) * s;
        d(k) = d(k) - z;
        d(k+1) = d(k+1) + z;
        od(k+1) = g * c * s + (c * c - s * s) * od(k+1);
        x = od(k+1);
        if (k > 1)
            od(k) = w;
        end
        if (k < (m-1))
            y = -s * od(k+2);
            od(k+2) = c * od(k+2);
        end
        U(:,k:(k+1)) = U(:,k:(k+1)) * [c s;-s c];
        
    end
    if ((abs(od(m)) < eps * (abs(d(m-1)) + abs(d(m)))) || (iter >= niter))
        if (iter >= niter)
            warning('ImplicitTriQR:instability',['ImplicitTriQR() did not converge for m = ' num2str(m)]);
        end
        m = m - 1;
        iter = 0;
    end
end
    D = diag(d);
    V = U;   
                
end

function dist = SafeDistance(a,b)
  abs_a = abs(a);
  abs_b = abs(b);
  if (abs_a > abs_b)
    dist = abs_a * sqrt(1.0 + (abs_b / abs_a)^2);
  else
    if (abs_b == 0)
      dist = 0;
    else
      dist = abs_b * sqrt(1.0 + (abs_a / abs_b)^2);
    end
  end
end

end
    

