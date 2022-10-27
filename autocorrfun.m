function [autocorr] = autocorrfun(F)
%[autocorr] = autocorrfun(F)
%   This function computes the autocorrelation of F. autocorr is the
%   autocorrelation function, which depends on the distance.
n=size(F,1);
m=size(F,2);
B=zeros(2^nextpow2(n)*2,2^nextpow2(m)*2);
for i=1:size(F,3)
    ft=[];
    ft=fft2(F(:,:,i),size(B,1),size(B,2));
    B= B+ft.*conj(ft);
end

B=ifftshift(ifft2(B));
B=B/B(n,m);
[autocorr] = AngAverage(B); % Angle average 
% function shown in the next code
end