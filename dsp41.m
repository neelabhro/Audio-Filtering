%Neelabhro Roy
%2016171

clc;
close all;
clear all;


%Question 1

%%%%%%%%% record the audio
Fs = 8000;
%Sampling Frequency
recObj = audiorecorder(Fs,16,1);
disp('Start speaking.')
recordblocking(recObj,5);
disp('End of Recording.');

% Play back the recording.
play(recObj);
% Store data in double-precision array.
hold on
hold on
hold on
hold on
myRecording = getaudiodata(recObj);
% Plot the waveform.
plot(myRecording);
title('Input Recording');
figure;

%Filtering
b = [0.25 0.5 0.25];
a=[1];
%n= 0:200;
%x = .1.^n;
% input is x
x = myRecording;

%d1 = fdesign.lowpass('N,Fc',10,1200,8000);
% We are performing Lowpass filtering at the Half Bandwidth
%designmethods(d1);

%f1 = design(d1, 'window');  
%fvtool(f);
y = filter(b,a,x);


%y = filter(b,a,x);
% filtered signal y
%[h w] = freqz(b,a,y);
plot(y);
title('Low Pass Filtered Output');
sound(y);
hold on
hold on
hold on
hold on
hold on
figure;
%plot(real(h),real(w));
freqz(y)
title('Magnitude Response of Low Pass filter');
figure;

%Question 2

%d2 = fdesign.highpass('N,Fc',10,1200,8000);
%('Fp,Fst,Ap,Ast',0.0001,0.1,1,60);
% We are doing Highpass Window filtering of the input signal, and as is
% evident from the resultant Sound, only the High frequencies above the
% Half Bandwidth are resolved and audible.
%designmethods(d2);
%f2 = design(d2, 'window');
b1 = [0.25 -0.5 0.25];
a1 = [1];

i = filter(b1,a1,x);
plot(i);
title('High Pass Filtered Output');
hold on
hold on
hold on
hold on
hold on
sound(i);
hold on
%[h1 w1] = freqz(b,a,i);
figure;
%plot(real(h1),real(w1));
freqz(i)
title('Magnitude Response of the High Pass filter');




figure;
n = 1024;
%N = length(y);
FFT = fft(x(1:1000),n);
stem(abs(FFT));
%Normalising the Fourier Transform
%Fn = Fs/2; 
% The Nyquist Frequency
%Freq = (( linspace(0,1,fix(N/2)+1)) .* Fn);
%Index = 1 : length(Freq);
%Assigning the indices
%stem(Freq./1000, abs(FFT(Index))*2);

title('Magnitude Spectrum of the Fourier Transform of the Audio signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
