%% Code For Bonus Assignment
% Registration Number D3200014
% Turkish ID = 
clear;
clc;

c1=6;
c2=5;
N = 2;


a1 = (0.89125)^2;
a2 = (0.1773)^2;

ourOmega = 2*tan(0.1*pi); 
ourOmega2 = 2*tan(0.15*pi);

OmegaC = (0.3+0.2)*pi/2;
tol = 10000;
while tol>1e-5
    b = 1/a1;
    c = b-1;
    d = log(c)/(2*N);
    OmegaCtemp = ourOmega/exp(d);
    b = 1/a2;
    c = b-1;
    N = log(c)/log(ourOmega2/OmegaCtemp)/2;
    tol = abs(OmegaC-OmegaCtemp);
    OmegaC = OmegaCtemp;
end
disp(ceil(N));
disp(OmegaC);

Nc=2*ceil(N);


%Finding Xc according to Number Of N
fun = @(x) a2-1/(1+((ourOmega2)/x)^Nc);  
x_true = fzero(fun,[0 1],optimset("Display","iter")); 

% finding poles of C(s)
a = 1;
b = (1/(x_true))^Nc;
VectorOfPolynomial = [b, zeros(1, Nc-1), a];
RootsOfPolynomial = roots(VectorOfPolynomial);

figure()
zplane([],RootsOfPolynomial)

RootsOfHalfZ=RootsOfPolynomial(real(RootsOfPolynomial)<0)
figure()
zplane([],RootsOfHalfZ)


K=RootsOfHalfZ(1)*RootsOfHalfZ(2)*RootsOfHalfZ(3)*RootsOfHalfZ(4)*RootsOfHalfZ(5)*RootsOfHalfZ(6)

b = K
a = poly(RootsOfHalfZ) 
[r, p, k] = residue(b,a)

poly([p(1) p(2)])
[num, den] = residue(r, exp(p), k)

[H,w]=freqz(num,den);
figure
plot(w/pi,20*log10(abs(H)));

syms s
syms z


Hs =b/expand((s-p(1))*(s-p(2)))*expand((s-p(3))*(s-p(4)))*expand((s-p(5))*(s-p(6)))
vpa(Hs,3)

s = 2*(1-z)/(1+z);

H1 = (s^2)
vpa(HS,3)

