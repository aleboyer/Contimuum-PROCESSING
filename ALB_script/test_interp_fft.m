%% test de l idee d interpolation via fft????

% create signal
x=linspace(-10*pi,10*pi,5000);
y=cos(2*x+5)-cos(10*x-6)+sin(3*x).^2;
plot(x,y)
legend('signal')
dx=mean(diff(x));
%put some nan in the middle
ynan=y;
ynan(400:405)=nan;ynan(1100:1200)=nan;
ynan(2000:2500)=nan;ynan(2000:2500)=nan;
hold on 
plot(x,ynan,'g','linewidth',2)
legend('signal nan')

yfilt=interp1(x(~isnan(ynan)),ynan(~isnan(ynan)),x);
hold on 
plot(x,yfilt,'r','linewidth',2)
legend('filt y')

ffty=fft(yfilt);
figure(2)
plot(abs(ffty))

% biggest nan pattern 500
bignan_freq=1/500;




