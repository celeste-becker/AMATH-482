%% PART 2 (Piano) Set Up - (lines 1-14)
clear; close all; clc

%set up
[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs; % total recorded time in seconds
t = (1:length(y))/Fs; % vector of time in seconds
S = y'; %Signal S equal to the amplitude  

L=tr_piano;
n=length(S); % n is the number of data points in signal  
k=(1/L)*[0:n/2-1 -n/2:-1]; %EVEN # of frequencys in Hz
ks=fftshift(k);

%% Spectrogram - (lines 15-57)
tslide=0:0.1:tr_piano; %Descrete time vector
Sgt_spec = zeros(length(tslide),n);

%defining gaussian Filter 
width = 0.00009;
centerF = 300;
filter = exp(-width*((ks - centerF).^2)); 

%gabor filter width
a = 150;

for j=1:length(tslide)
    
    %applying filter
    g=exp(-a*(t-tslide(j)).^2); 
    Sg=g.*S; 
    Sgt=fft(Sg);%fourier transforming to the freuqency domain
    
    %Gaussian filter function --> use non-shifted k because the sgt uses
    %non-shifted k
    filter = exp(-width*((k - centerF).^2)); 
    
    Sgt = Sgt .* filter;
    
    %storing all the filtered signal at that time in the spectogram vector
    Sgt_spec(j,:) = fftshift(abs(Sgt)); 

end


%plot spectogram
pcolor(tslide,ks,Sgt_spec.'), 
shading interp 
title('Piano Spectogram');
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
set(gca,'Fontsize',16) 
colormap(hot)
colorbar
ylim([240 340]);
p8 = audioplayer(y,Fs); playblocking(p8);
