function [E,varargout]=KinEnSpectrum(vx,vy,x,y)
% [q,Eq]=KinEnSpectrum(vx,vy,Mx,My)
% INPUT:
% --> vx, vy (matrices containing x and y coordinate of the velocity.
% They can be 3d matrices with all the frames in the 3rd
% dimension of the matrices)
% --> x, y position matrices
% OUTPUT:
% E: 1d power spectrum of the kinetic energy
% q (optional)



Nframes=size(vx,3);
Nx=size(x,2); Ny=size(y,1);

Fsx=1/(x(1,2)-x(1,1)); % sampling frequency along x direction
Fsy=1/(y(2,1)-y(1,1)); % sampling frequency along y direction
Nf=2.^(nextpow2([Ny,Nx])); % size of the output of FFT
qx=linspace(0,Fsx/2,Nf(2)/2)*2*pi; % wave-number along x
qy=linspace(0,Fsy/2,Nf(1)/2)*2*pi; % wave-number along y
[Mqx,Mqy]=meshgrid([-qx(end:-1:1),qx],[-qy(end:-1:1),qy]);

% compute the 2D kinetic energy spectrum for each frame
for iframe=1:Nframes
    V(:,:,1)=vx(:,:,iframe); V(:,:,2)=vy(:,:,iframe);
    for ic=1:2 % compute the FFT
        ftV(:,:,ic)=fftshift(fft2(V(:,:,ic),Nf(1),Nf(2)))/(Nx*Ny);
    end
    ME2d(:,:,iframe)=(abs(ftV(:,:,1)).^2+abs(ftV(:,:,2)).^2)/2;
end
Lx=abs(x(end)-x(1));
Ly=abs(y(end)-y(1)); %system size (Lx*Ly)
E2d=mean(ME2d,3)*Lx*Ly; % Temporal average
Mq=sqrt(Mqx.^2+Mqy.^2);
E=AngularAverage(E2d.*Mq); %Angular average after multiplying by Jacobian
if nargout>1
    q=AngularAverage(Mq);
    varargout{1}=q;
end

    function [avgF,varargout]=AngularAverage(F,varargin)
        %   [avgF,varargout]=AngularAverage(F,varargin)
        % INPUT: F (2d matrix to be angle averaged)
        % optional input: x and y coordinates
        % OUTPUT: avgF (vector containing the angle average of F
        % optional output: vector with the distance from the center

        Nx=size(F,2); Ny=size(F,1);
        imp_coord=0;
        if nargin>1
            Mposx=varargin{1}; Mposy=varargin{2};
            imp_coord=1;
        end

        if ~imp_coord
            [Mposx,Mposy]=meshgrid(0:Nx-1,0:Ny-1);
        end

        if mod(Nx,2)==0
            centerx=(Mposx(1,floor((Nx+1)/2))+Mposx(1,ceil((Nx+1)/2),1))/2;
        else
            centerx=Mposx(1,(Nx+1)/2);
        end

        if mod(Ny,2)==0
            centery=(Mposy(floor((Ny+1)/2),1)+Mposy(ceil((Ny+1)/2),1))/2;
        else
            centery=Mposy((Ny+1)/2,1);
        end
        Mdistx=Mposx-centerx; Mdisty=Mposy-centery;
        Mr=sqrt(Mdistx.^2+Mdisty.^2);
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
        bins=discretize(Mr,r);

        for ibin=1:numel(r)-1
            avgF(ibin)=mean(F(bins==ibin));
            if nargout>1
                R(ibin)=mean(Mr(bins==ibin));
            end
        end
        if nargout>1
            varargout{1}=R;
        end
    end

end