N = 512;
x = linspace(-1.5,1.5,N)

sigh = 0.24;
%H = exp(-(x).^2/(2*sigh^2));
Ampdom = [.61 1 .61];
mudom = [-sigh 0 sigh];

f1 = 1;
f2 = 4;

sig1 = 1/(pi*f1);
sig2 = 1/(pi*f2);

figure
for i = 1:length(mudom)
    
    mu = mudom(i);
    amp = Ampdom(i);
    
    subplot(1,2,1)
    hold on
    y = amp*exp(-(x-mu).^2/(2*sig1^2));
    plot(x,y,'b')
    axis off
    
    subplot(1,2,2)
    hold on
    y = amp*exp(-(x-mu).^2/(2*sig2^2));
    plot(x,y,'b')
    axis off
 
end

sigp1 = sqrt(sig1^2 + sigh^2);
subplot(1,2,1)
hold on
y = 1.05*exp(-(x).^2/(2*sigp1^2));
plot(x,y,'g')
axis off

sigp2 = sqrt(sig2^2 + sigh^2);
subplot(1,2,2)
hold on
y = 1.05*exp(-(x).^2/(2*sigp2^2));
plot(x,y,'g')
axis off

%%
N = 512;
x = linspace(-8,8,N);

sigh = 0.85;
%H = exp(-(x).^2/(2*sigh^2));
Ampdom = [.61 1 .61];
mudom = [-sigh 0 sigh];

f1 = 1;
f2 = 4;

sig1 = f1/2;
sig2 = f2/2;

figure
for i = 1:length(mudom)
    
    mu = mudom(i);
    amp = Ampdom(i);
    
    subplot(1,2,1)
    hold on
    y = amp*exp(-(x-mu).^2/(2*sig1^2));
    plot(x,y,'b')
    axis off
    
    subplot(1,2,2)
    hold on
    y = amp*exp(-(x-mu).^2/(2*sig2^2));
    plot(x,y,'b')
    axis off
 
end

sigp1 = sqrt(sig1^2 + sigh^2)
subplot(1,2,1)
hold on
y = 1.05*exp(-(x).^2/(2*sigp1^2));
plot(x,y,'g')
axis off
%xlim([-4 4])

sigp2 = sqrt(sig2^2 + sigh^2)
subplot(1,2,2)
hold on
y = 1.05*exp(-(x).^2/(2*sigp2^2));
plot(x,y,'g')
axis off

