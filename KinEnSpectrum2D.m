function [E2d,Mqx,Mqy]=KinEnSpectrum2D(vx,vy,x,y)
% [E2d,Mqx,Mqy]=KinEnSpectrum(vx,vy,Mx,My)
% INPUT: 
% --> vx, vy (matrices containing x and y coordinate of the velocity. 
% They can be 3d matrices with all the frames in the 3rd 
% dimension of the matrices)
% --> x, y position matrices
% OUTPUT:
% E2d: 2d power spectrum of the kinetic energy
% Mqx, Mqy wavenumbers


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
end