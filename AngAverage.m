function [avgF,stdF,varargout]=AngAverage(F,varargin)
% This function computes the angle average of a 2D matrix
%   [avgF,varargout]=AngularAverage(F,varargin)
% INPUT: F (2d matrix to be angle averaged)
% optional input: x and y coordinates and angles between which
% we can perform the average, default: 0 --> 2pi
% OUTPUT: avgF and stdF (vectors containing the mean and the std of F)
% optional output: vector with the distance from the center

Nx=size(F,2); Ny=size(F,1); Ntime=size(F,3);
lim_angles=0;
imp_coord=0;

if nargin>1
    for ivar=1:2:numel(varargin)
        if length(varargin{ivar})==1
            angle1=varargin{ivar}; angle2=varargin{ivar+1};
            lim_angles=1;
        else
            Mposx=varargin{ivar}; Mposy=varargin{ivar+1};
            imp_coord=1;
        end
    end
end

if ~imp_coord
    if mod(Nx,2)==0
%         Fbig=[F,zeros(Ny,1)];
%         F=Fbig;
        F(:,end+1,:)=0;
        Nx=size(F,2);
    end
    if mod(Ny,2)==0
%         Fbig=[F;zeros(1,Nx)];
%         F=Fbig;
        F(end+1,:,:)=0;
        Ny=size(F,1);
    end
    x=linspace(-size(F,2)/2,size(F,2)/2,Nx); y=linspace(size(F,1)/2,-size(F,1)/2,Ny); 
    [Mposx,Mposy]=meshgrid(x,y);
end


if lim_angles
    angles=atan2(Mposy,Mposx);
    if double(angle1<0)+double(angle2<0)~=1
        angles=angles+2*pi*(angles<0);
    end
end

Mr=sqrt(Mposx.^2+Mposy.^2);
if ~lim_angles
    if Nx>Ny
        if mod(Nx,2)==0
            r=Mr(round(Ny/2),round(Nx/2)+1:end);
        else
            r=Mr(round(Ny/2),round(Nx/2):end-1);
        end
    else
        if mod(Ny,2)==0
            r=Mr(round(Ny/2)+1:end,round(Nx/2))';
        else
            r=Mr(round(Ny/2):end-1,round(Nx/2))';
        end
    end
else
    if Nx>Ny
        rmax=max(max(Mposx));
    else
        rmax=max(max(Mposy));
    end
    Nbins=round(max([Nx,Ny])/2);
    [N,r]=histcounts(Mr(Mr<=rmax&angles>=angle1&angles<=angle2),Nbins);
    while min(N)==0
        Nbins=Nbins-1;
        [N,r]=histcounts(Mr(Mr<=rmax&angles>=angle1&angles<=angle2),Nbins);
    end
end
bins=discretize(Mr,r);

for ibin=1:numel(r)-1
    if ~lim_angles
        if Ntime>1
            [i1,i2]=find(bins==ibin);
            i1frames=repmat(i1',1,Ntime);
            i2frames=repmat(i2',1,Ntime);
            i3frames=reshape(repmat(1:Ntime,numel(i1),1)...
            ,[1,numel(i2frames)]);
            itake=sub2ind([Ny,Nx,Ntime],i1frames,i2frames,...
            i3frames);
            stdEll(ibin)=std(F(itake));
            avgF(ibin)=mean(F(itake));
        else
            avgF(ibin)=mean(F(bins==ibin));
            stdEll(ibin)=std(F(bins==ibin));
        end
    else
        if Ntime>1
            [i1,i2]=find(bins==ibin&angles>=angle1&...
            angles<=angle2);
            i1frames=repmat(i1',1,Ntime);
            i2frames=repmat(i2',1,Ntime);
            i3frames=reshape(repmat(1:Ntime,numel(i1),1),...
            [1,numel(i2frames)]);
            itake=sub2ind(size(F),i1frames,i2frames,i3frames);
            stdEll(ibin)=std(F(itake));
            avgF(ibin)=mean(F(itake));
        else
            stdEll(ibin)=std(F(bins==ibin&angles>=angle1...
            &angles<=angle2));
            avgF(ibin)=mean(F(bins==ibin&angles>=angle1...
            &angles<=angle2));
        end
    end
end
varargout{1}=(r(1:end-1)+r(2:end))/2;
end