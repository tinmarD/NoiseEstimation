

t   = 0:1/2000:2-1/2000;
x   = 1+chirp(t-2,4,1/2,6,'quadratic',100,'convex').*exp(-4*(t-1).^2);
%- Calcul the mean of the signal
xMean = mean(x);
%- Compute the hilbert transform of the signal minus the mean
env = abs(hilbert(x-xMean));
%- Add the mean to the envelope
yUp  = env+xMean;
yLow = xMean-env;
%- Median of the envelope gives an estimate of the width of the signal
s_envMed  = median(abs([yUp,yLow]));
disp(['Epaisseur du signal : ',num2str(s_envMed)]);


figure; hold on;
plot(t,x);
plot(t,yUp,'r','linewidth',2)
plot(t,yLow,'g','linewidth',2);


% t = 0:1e-4:1;
% x = [1+cos(2*pi*50*t)].*cos(2*pi*1000*t);
% 
% y = hilbert(x);
% env = abs(y);
% 
% figure;
% plot(t,x)
% hold on
% plot(t,[-1;1]*env,'r','LineWidth',2)
% xlim([0 0.1])
% xlabel('Seconds')
