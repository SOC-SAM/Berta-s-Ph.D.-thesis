function [filtdirector]=FilterDirector(director, filtersize)
%[filtdirector]=FilterDirector(director, filtersize) 
%   This function applies a disk filter (averages the values over a disk
%   of radius filtersize) to an nematic orientational field (director)
%   and gives filterdirector

    filt=fspecial('disk',filtersize/2); % Create the filter
    Qxx=cos(2*director)/2; Qxy=sin(2*director)/2; 
    % Apply the filter
    filtQxx=imfilter(Qxx,filt,'conv'); 
    filtQxy=imfilter(Qxy,filt,'conv');
    % Compute S and \theta
    filtSS=2*sqrt(filtQxx.^2+filtQxy.^2);
    v1=ones(size(filtSS));
    v2=(-filtQxx+filtSS/2);
    filtdirector=atan2(v2./filtQxy,v1); 
    filtdirector=filtdirector+pi*(filtdirector<0); % to make sure the angles go from 0 to pi
end