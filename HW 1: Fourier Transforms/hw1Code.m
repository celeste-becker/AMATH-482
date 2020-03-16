clear; close all; clc;

%% Load Data - (lines 3-5)
load Testdata

%% Set up - (lines 6-13)
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

%% Average in Frequency Domain - (lines 14-24)
avg = zeros(n,n,n);
for row = 1:20
    Un = reshape(Undata(row,:),n,n,n);
    
    Unf = fftn(Un);
    
    avg = avg + Unf;
end
avg = abs(avg./20);

%% Find Center Frequency From Averaged Data - (lines 25-32)
[maxValue, indexOfMax] = max(avg,[],'all', 'linear'); 

%center frequencys in Kx, Ky and Kz
Kx0 = Kx(indexOfMax);
Ky0 = Ky(indexOfMax);
Kz0 = Kz(indexOfMax);

%% Gaussian Filtering - (lines 33-63)
%Tau = the bandwidth of the filter 
t = 1;

%Gaussian filter function
filter = exp(-t*((Kx-Kx0).^2 + (Ky-Ky0).^2 + (Ky - Ky0).^2));

%Marble Position vectors
x0 = zeros(1, 20);
y0 = zeros(1, 20);
z0 = zeros(1, 20);

for row = 1:20
    
    Un = reshape(Undata(row,:),n,n,n);
    
    Unf = fftn(Un); 
    
    UfilterFS = Unf.*filter;
    
    UfilterSS = abs(ifftn(UfilterFS));  
    
    [maxPos, I] = max(UfilterSS(:));
    [i1, i2, i3] = ind2sub([n,n,n], I);
    
    x0(1, row) = X(i1, i2, i3);
    y0(1, row) = Y(i1, i2, i3);
    z0(1, row) = Z(i1, i2, i3); 
     
end

%% Plot Marble Trajectory and Last Data Point - (lines 64-78)
plot3(x0,y0,z0, 'r:', 'Linewidth', [3]);
hold on
plot3(x0(20), y0(20), z0(20),'b*', 'Linewidth', [10]);
title('Trajectory of Marble');
xlabel('X axis');
ylabel('Y axis');
zlabel('z axis');