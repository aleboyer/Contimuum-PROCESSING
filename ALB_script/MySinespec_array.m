function [f,P,err]=MySinespec_array(data,samplefreq,K);
%   [f,P]=MySineSpec(data,samplefreq,K) takes an input velocity 
%   vector 'data' (complex), sample frequency in samples per day, K 
%   multitapers (generally 10), and returns a sine multitaper 
%   spectrum P with frequency vector f.   
%   
%   when plotting, loglog(f,P,-f,P) to see both sides.

%   See also MySineTaper (real)  MW 8/03
%   Change to make it work with 2d array ALB 01/13/2016
err=0;
if ~isempty(find(isnan(data(:))))
    err=1;
end
[Z,N]=size(data);
if length(data)~=N;
    disp(num2str(Z))
    disp(num2str(N))
    disp('time < Z is it normal?')
end
nfft=N;
u=real(data);
v=imag(data);
u=detrend(u','constant')';
v=detrend(v','constant')';
rex=complex(u,v);
A=ones(N,Z,K);
for k=1:K;
	A(:,:,k)=((2/(N+1))^(1/2)*sin(((k)*pi*(1:N))/(N+1)))'*ones(1,Z);
	A(:,:,k)=A(:,:,k).*rex';  %this gives us h*X in the formula.  
end
fta=fft(A,[],1);
msftA=(1/samplefreq)*(abs(fta)).^2;
sumk=sum(msftA,3)./K;  

dz=1/samplefreq;
dk=(1/N)/dz;
%This is not quite right if N is odd...  fix.  mha 6/25/99
if rem(N,2)==0
k=-N/2*dk:dk:N/2*dk-dk;
else
	kp=dk:dk:dk*floor(N/2);
	k=[fliplr(-kp) 0 kp];
end

f=k;
Pa=sumk';
mid=ceil(length(k)/2);
P1=Pa(:,mid:length(Pa));
P2=Pa(:,2:mid);
P=[P1 P2];
