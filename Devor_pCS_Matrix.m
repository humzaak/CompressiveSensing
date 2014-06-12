function Phi = Devor_pCS_Matrix(n , k , q)




% Phi = Devor_pCS_Matrix(n , k)
%
% This function constructs a (q^2 * n) binary matrix (prior to column 
% normalization) which has a coherence less than 1/(k-1). This is a 
% sufficient condition for RIP of order "k" and also perfect reconstrcution 
% of (k/2)-sparse vectors. The design is based on the polynomials in Galois 
% fields (Devore scheme) and the details can be found in:
%
% [*] R. DeVorea"Deterministic constructions of compressed sensing matrices,"
%     Journal of Complexity, vol. 23, no. 4-6, Aug.-Dec. 2007, pp. 918-925.
%
% "k" and "n" are the standard notations of the compressed sensing parameters. 
%
% "q" is the size of the Galois Field and should be a prime power.
%
% The function is written by Arash Amini


% check if q is a pime power
q_factors       = factor(q);
if q_factors(1) == q_factors(end)
    p           = q_factors(1);
    power       = length(q_factors);
else
    error('Error:   "q" must be a prime power')
end


% finding polynomial degree
r               = floor(q / (k-1));
if r * (k-1) == q
    r 	    = r - 1;
end

% compatibilty check
if n > q ^ (r+1)
    disp(['Warning:  using the input "k" and "q", the maximum reachable "n" is ' , num2str(q ^ (r+1))])
elseif n < q ^ (r+1)
    disp(['!!!! Note: The same number of samples can be used (the same "k") for n <= ' , num2str(q ^ (r+1)) , ' !!!!'])
end


% setting the field and its primitive root
field           = gftuple((-1 : q - 2).' , power , p); %gftuple only produces a simplified polynomial-
alpha           = 2;        % the primitive root




% constrcuting the matrix by using different polynomials for different columns
Phi             = zeros(q ^ 2  ,  n);
for Poly_ind = 1 : n
    PolyCoeff   = de2bi(Poly_ind - 1  ,  r + 1  ,  q , 'left-msb') - 1; %de2bi(d,n,p,'flg')
    
    PolyVal     = -ones(1 , q);
    for coef_ind = 1 : r+1
        PolyVal = gfmul(PolyVal , [-1 : q - 2] , field);
        PolyVal = gfadd(PolyVal , repmat(PolyCoeff(coef_ind) , 1 , q) , field); %add normally and then take modulo
    end

    PolyVal(PolyVal == -inf)    = -1;
    Phi(2 + PolyVal + q * [0 : q - 1]  ,  Poly_ind) = 1;
end


% normalization
assignin('base','Phi',Phi);
% Phi             = Phi / sqrt(q);