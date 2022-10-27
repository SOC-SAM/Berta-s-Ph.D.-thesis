function [WindingNumber] = winding(theta,center,WindowSize,varargin)
%[WindingNumber] = winding(theta,center,WindowSize
% This function computes the WindindNumber of theta nematic director field
% around center [row, col] in MATRIX INDEXES NOTATION;
% theta by default is referenced to pixel indices coordinate system.
%(0,0) --> Upper left pixel.
% If theta must be considered in cartesian reference system call winding:
% [WindingNumber] = winding(theta,center,WindowSize,'cartesian')
Nx=size(theta,2); Ny=size(theta,1);
[ ind ] = index_circle( [center(2),center(1)],WindowSize,Nx,Ny, 0, 2*pi);


for i=1:size(ind,1)-1
    dif(i)=theta(ind(i+1,1),ind(i+1,2))-theta(ind(i,1),ind(i,2));
    if abs(dif(i))>pi/2
        dif(i)=theta(ind(i+1,1),ind(i+1,2))-...
            theta(ind(i,1),ind(i,2))+(sign(dif(i)))*pi;
    end
end
WindingNumber=sum(dif)/(2*pi);
if isempty(varargin)==0
    if strcmp(varargin{1},'cartesian')
        WindingNumber=sum(dif)/(2*pi)*-1;
    end
end


    function [ ind ] = index_circle( center,size_cir,Nx,Ny, initial_ang, final_ang)
        %INDEXES IN A CIRCLE
        %   This function computes the indexes of a matrix in a circle. It doesn't
        %   work for large size_cir.
        % [ ind ] = index_circle( center,size_cir,Nx,Ny, initial_ang, final_ang)

        comp_ind=Inf;
        comp_ind2=Inf;
        l1=1;
        l2=1;
        dphi=(2*pi)./(2*pi*size_cir*1.4);
        angle=(initial_ang:dphi:final_ang);
        for i=1:size(angle,2)
            kk(1)=round(cos(angle(i))*size_cir+center(1));
            kk(2)=round(sin(angle(i))*size_cir+center(2));
            logi1=kk~=comp_ind;
            if max(logi1)
                comp_ind=kk;
                ind(l1,:)=kk;
                l1=l1+1;
            end
        end
    end
end