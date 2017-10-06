%%NLMS CALGORITHM FOR THE SOURSE CODE
%%
function [A,E,Y] = nlms(x,d,beta,nord,a0)
X=convm(x,nord);

[M,N]=size(X);
if nargin < 5, a0 = zeros(1,N); end
%initialization
a0=a0(:).';
Y(1)=a0*X(1,:).';  
E(1)=d(1) - a0*X(1,:).';
DEN=X(1,:)*X(1,:)'+0.0001;

A(1,:) = a0 + beta/DEN*E(1)*conj(X(1,:));
if M>1
    for k=2:M-nord+1;
        Y(k)=A(k-1,:)*X(k,:).';%output equation
        E(k) = d(k) - A(k-1,:)*X(k,:).';%error signal
        DEN=X(k,:)*X(k,:)'+0.0001;%normalizing the input signal
        A(k,:)=A(k-1,:)+ beta/DEN*E(k)*conj(X(k,:));%update equation
        
    end;
end;
