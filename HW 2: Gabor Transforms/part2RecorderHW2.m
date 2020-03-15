%% Part 2 (Recorder) - (lines 1-20)
clear; close all; clc
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs; % total recorded time in seconds
t = (1:length(y))/Fs; % vector of times
S = y'; %sets the signal equal to y (amplitude) of the filter 

%setting frequency domain stuff
L=tr_rec;
n=length(S); % n is the number of data points in signal  
k=(1/L)*[0:n/2-1 -n/2:-1]; %EVEN # of frequencys **** DEVIDE BY L OR NAH???
% ^ by changing the scaling factor for k to (1/L) this gives the freuqencys
% in Hz
ks=fftshift(k);

%gaussian Filter --> used to filter out overtones 
width = 0.00005;
centerF = 900;
filter = exp(-width*((ks - centerF).^2)); 

%% Spectrogram - (lines 21-60)

tslide=0:0.1:tr_rec; % this makes the time vector in steps
Sgt_spec = zeros(length(tslide),n); % defines the spectogram freuqency vector 
%^ sgt_spec(time in increments of tslide, length of signal)

%gabor filter width
a = 150;

for j=1:length(tslide)
    
    %applying filter
    g=exp(-a*(t-tslide(j)).^2); 
    Sg=g.*S; 
    %fourier transforming to the freuqency domain
    Sgt=fft(Sg); 
    
    %Gaussian filter function --> use non-shifted k because the sgt uses
    %non-shifted k
    filter = exp(-width*((k - centerF).^2)); 
    
    Sgt = Sgt .* filter;
    
    %storing all the filtered signal at that time in the spectogram vector
    Sgt_spec(j,:) = fftshift(abs(Sgt)); 
end

%plot spectogram
figure(2)
pcolor(tslide,ks,Sgt_spec.'), 
shading interp 
title('Recorder Spectogram');
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
set(gca,'Fontsize',16) 
colormap(hot)
colorbar
ylim([775 1050]);

p8 = audioplayer(y,Fs); playblocking(p8);