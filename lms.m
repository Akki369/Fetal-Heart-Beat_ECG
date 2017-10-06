%%LMS CODE ALGORITHM FOR THE SOURSE CODE
%%
function [A,E,Y] = lms(x,d,mu,nord)
%lms function
X=convm(x,nord);
[M,N]=size(X);
if nargin < 5, a0 = zeros(1,N); end
a0=a0(:).';

Y(1)=a0*X(1,:).';
E(1)=d(1) - Y(1);

A(1,:) = a0 + mu*E(1)*conj(X(1,:));
if M>1
    for k=2:M-nord+1;
        Y(k,:)=A(k-1,:)*X(k,:).';%ouput equation
        E(k,:) = d(k) - Y(k,:);%error signal
        A(k,:)=A(k-1,:)+mu*E(k)*conj(X(k,:));%update equation
    end;
end;
%%
