clc
clear all
close all
%%
N=96;
L=5;
fs=N/L;
t = linspace(0,L,N);
%y = chirp(t,5,L,15);
y = sin(2*pi*10*t);
spectrogram(y,64,32,32,fs)

%%
fft_size=64;
F=zeros(fft_size,fft_size);
w = hann(fft_size);

for k=1:fft_size
    for j=1:fft_size
        
        F_h(k,j)=w(j)*(exp(1i*j*k*2*pi/fft_size)/sqrt(fft_size));
    end
    
end

A=zeros((N/32-1)*fft_size,N);

x_shift=32;
y_shift=fft_size;

for k=1:(N/32-1)
        current_x=x_shift*k;
        current_y=y_shift*k;
        A(y_shift*k-fft_size+1:y_shift*k,x_shift*k-fft_size/2+1:x_shift*k+fft_size/2)=F_h;
    
end

figure (2)
surf(real(A))
xlim([0 N])
ylim([0 (N/32-1)*fft_size])
set ( gca, 'xdir', 'reverse' )
view(2)
x0=10;
y0=10;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

%%
F_cal=A*y';
F_mat=reshape(F_cal,64,N/32-1);
F_mat_single=F_mat(1:32,:)*2/sqrt(fft_size);
%%
fm=linspace(0,fs/2,fft_size/2);
tm=linspace(0,L,N/32-1);

[TM,FM] = meshgrid(tm,fm);

figure (3)
surf(TM,FM,(10*log10(F_mat_single.*conj(F_mat_single))))
xlim([0 L])
ylim([0 fs/2])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
view(2)
grid off
colorbar

%% HW 1 Extra

clc
clear all
close all
%%

recObj = audiorecorder
disp('Start speaking.')
recordblocking(recObj, 1);
disp('End of Recording.');

%%
play(recObj);
y = getaudiodata(recObj);
%%
t=linspace(0,1,8000);
plot(y);

%%
y_new=y(2049+768:7168);
fs=8000;

t_new=linspace(0,(7168-2049-768)/fs,7168-2049-768+1);
N=size(y_new,1);

spectrogram(y_new,256,128,128,fs)

%%
fft_size=256;
F=zeros(fft_size,fft_size);
w = hann(fft_size);

for k=1:fft_size
    for j=1:fft_size
        
        F_h(k,j)=w(j)*(exp(1i*j*k*2*pi/fft_size)/sqrt(fft_size));
    end
    
end

A=zeros((N/128-1)*fft_size,N);

x_shift=128;
y_shift=fft_size;

for k=1:(N/128-1)
        current_x=x_shift*k;
        current_y=y_shift*k;
        A(y_shift*k-fft_size+1:y_shift*k,x_shift*k-fft_size/2+1:x_shift*k+fft_size/2)=F_h;
    
end

figure (2)
surf(real(A))
xlim([0 N])
ylim([0 (N/128-1)*fft_size])
set ( gca, 'xdir', 'reverse' )
view(2)
x0=10;
y0=10;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

%%
F_cal=A*y_new;
F_mat=reshape(F_cal,fft_size,N/128-1);
F_mat_single=F_mat(1:128,:)*2/sqrt(fft_size);
%%
L=t_new(end);
fm=linspace(0,fs/2,fft_size/2);
tm=linspace(0,L,N/128-1);

[TM,FM] = meshgrid(tm,fm);

figure (4)
surf(TM,FM,(10*log10(F_mat_single.*conj(F_mat_single))),'edgecolor','none')
xlim([0 L])
ylim([0 fs/2])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
view(2)
grid off
colorbar
set ( gca, 'xdir', 'reverse' )

