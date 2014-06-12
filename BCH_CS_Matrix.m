function Phi = BCH_CS_Matrix(n , k)

n  = 16
k = 4




% Phi = BCH_CS_Matrix(n , k)
%
% This function constructs an RIP fulfiling bipolar matrix which is based 
% on the BCH-based design of the following paper:
%
% [*] A. Amini and F. Marvasti, ?Deterministic Construction of Binary, 
%     Bipolar and Ternary Compressed Sensing Matrices," IEEE Trans. on 
%     Inform. Theory, vol. 57, no. 4, pp. 2360-2370, April 2011.
%
% The generated matrix is an m*n column-normalized bipolar matrix with a
% coherence value less than 1/(k-1) which is a sufficient condition for 
% establishing the RIP of order "k" and also perfect reconstrcution of 
% (k/2)-sparse vectors. In fact, "k", "n" and "m" are the standard notations 
% in compressed sensing.
%
% The function is written by Arash Amini



% parameters
RootPower       = 1;        % choosing the primitive root of the field by 
                            % setting a power for the default primitive root (2) in MATLAB

% determining \tilde{m} for the requested n
i_alg           = ceil(log2(k));
tildem          = i_alg - 1;
hx              = [1];
 
while 2 ^ (length(hx) - 2) < n
    tildem      = tildem + 1;
    [hx , gx , hxtext]  = PolyBCH(tildem , i_alg , RootPower);
end
m               = 2 ^ tildem  -1;

if 2 ^ (length(hx) - 2) > n
    disp(['!!!! Note: The same number of samples can be used (the same "k") for n <= ' , num2str(2 ^ (length(hx) - 2)) , ' !!!!'])
end


% modifing the code for the even parity
gx              = mod([gx  0]  + [0  gx]  ,  2);
gx              = double(gx);


% generating the matrix
Phi             = zeros(m  ,  n);
for aa =0 : n-1 
%     aa
    uncoded     = de2bi(aa  ,  'left-msb'  ,  length(hx) - 2);
    Phi(:  ,  aa + 1) = mod(conv(uncoded , gx)  ,  2).';
end


% Bipolar representation and normalization
% Phi(Phi == 0)   = -1;
%Phi             = Phi / sqrt(m);







function [hx , gx , hxtext] = PolyBCH(tildem , i , RootPower)

% [hx , gx , hxtext] = PolyBCH(tildem , i , RootPower)
%
% This function finds the binary parity check ("hx") and code generating 
% ("gx") polynomials of the BCH code used for CS purposes.
%
% The roots of these binary polynomials belong to GF(2^tildem); moreover, 
% the roots of h(x) are restricted to lie in the set 
% {1, alpha, alpha^2, ..., alpha^( 2^(tildem-1) + 2^(tildem - i) )}
% where "alpha" is one of the primitive roots of the field.
%
% The input "RootPower" is an optional input which specifies the primitive
% root. The default value for the root is '2'; however, if "RootPower" is
% given, this root is 2^RootPower.
%
% Note:     "RootPower" should be prime relative to 2^tildem - 1; when this
%           condition is violated, an error mesage will appear.
%
% Remark:   The outputs "hx" and "gx" are binary row vectors which show the 
%           polynomial coefficients with the descending order with respect 
%           to the powers of 'x'; i.e., the last elements are constant 
%           terms of the polynomial. In addition, "hxtext" is a string
%           which shows the parity check polynomial.
%
% The function is written by Arash Amini

% setting the primitive root
if nargin == 2
    alpha       = gf(2 , tildem);
else
    if gcd(RootPower , 2^tildem - 1) == 1
        alpha   = gf(2 , tildem) ^ RootPower;
    else
        disp(['Error:    The input ''RootPower'' should be prime to 2^tildem - 1 which is ' , num2str(2^tildem - 1) , ' in this setting'])
        return
    end
end

% finding the root powers of the parity check polynomial ( h(x) )
Hbin            = binHseq(tildem , i);
Hseq            = bi2de(Hbin , 'left-msb');

% evaluating the parity check polynomial by its roots
Hpoly           = gf(1 , tildem);
for bb = 1 : length(Hseq)
    Hpoly       = conv(Hpoly  , [alpha^0  alpha^(Hseq(bb))]);
end
hx              = uint8(Hpoly.x);
deg_h           = length(hx) - 1;

% producing the string that represents the parity check polynomial
hxtext          = ['     h(x) = 1'];
for bb = deg_h : -1 : 1
    if hx(bb) == 1
        hxtext  = [hxtext , ' + x^' , num2str(length(Hseq) - bb + 1)];
    end
end

% finding g(x) by dividing x^(2^tildem - 1) + 1 by h(x)
BigPoly         = uint8([1    zeros(1  ,  2^tildem - 2)     1]);
gx              = uint8(zeros(1  ,  2^tildem - deg_h));
for bb = 1 : length(gx)
    if BigPoly(bb) == 1
        BigPoly(bb : bb + deg_h)    = mod(BigPoly(bb : bb + deg_h) + hx  ,  2);
        gx(bb)  = 1;
    end
end









function Hbin = binHseq(tildem , i)

% Hbin = binHseq(tildem , i)
%
% This function finds all binary sequences of length "tildem" such that 1s
% are circularly spaced with at least $i$ zeros. These binary sequences are
% put together as rows of a matrix ("H")
%
% The function is written by Arash Amini

Allseq          = de2bi([0 : 2^tildem - 1]  ,  'left-msb');
isOK            = zeros(1 , 2^tildem);
AllInd          = 1 : tildem;

for aa = 1 : 2^tildem
    seq         = Allseq(aa , :);
    loc1        = AllInd(seq == 1);
    
    
    if length(loc1) > 1         
        diffs   = loc1(2 : end) - loc1(1 : end-1);          % spacing between the 1s in the middle
        diffs   = [loc1(1) + tildem - loc1(end)    diffs];  % including the circular spacing between 
                                                            % the first and the last 1
        if min(diffs) > i       % whether there exists at least "i" zeros in between
            isOK(aa)    = 1;
        end
        
    else
        isOK(aa)= 1;
    end
end

Hbin            = Allseq(isOK == 1  ,  :);