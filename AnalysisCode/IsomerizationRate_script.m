pF0 %Set paths
%%

[pRodDaylight IsomRateDaylight elevmu pupDiDom] = getIsomerizationRate_daylight(1);

[pRodsCRT IsomRateCRT LumDom pupDiDom] = getIsomerizationRate_LCD(1);

IsomDom = logspace(log10(10),log10(10^7),100);
param = [.5863 .0722];
pRodsModel = exp(-param(2)*(IsomDom/1000).^param(1));


%%

for i = 1:length(pupDiDom)
    legstr{i} = num2str(pupDiDom(i));
end

IsomtickDom = logspace(log10(10^0),log10(10^7),8);
LumtickDom = logspace(log10(10^0),log10(10^3),4);
ElevtickDom = -30:30:90;


%CRT %%%%%%
figure,
subplot(2,3,1)
loglog(LumDom,IsomRateCRT)
ylabel('R*/rod/sec')
xlabel('Luminance (cd/m^2)')
set(gca,'XTick',LumtickDom)
set(gca,'YTick',IsomtickDom)
ylim([IsomtickDom(1) IsomtickDom(end)])
xlim([LumtickDom(1) LumtickDom(end)])

subplot(2,3,4)
semilogx(LumDom,pRodsCRT)
ylabel('rods/(rods+cones)')
xlabel('Luminance (cd/m^2)')
ylim([0 1])
xlim([LumtickDom(1) LumtickDom(end)])
set(gca,'XTick',LumtickDom)
hold on 
plot([50 50],[0 1],'--b')

%Solar%%%%%%%%
subplot(2,3,2)
semilogy(elevmu,IsomRateDaylight)
ylabel('R*/rod/sec')
xlabel('solar elevation (deg)')
xlim([-30 90])
set(gca,'YTick',IsomtickDom)
ylim([IsomtickDom(1) IsomtickDom(end)])
set(gca,'XTick',ElevtickDom)
xlim([ElevtickDom(1) ElevtickDom(end)])

subplot(2,3,5)
plot(elevmu,pRodDaylight)
ylabel('rods/(rods+cones)')
xlabel('solar elevation (deg)')
ylim([0 1])
set(gca,'XTick',ElevtickDom)
xlim([ElevtickDom(1) ElevtickDom(end)])
%set(gca,'YTick',[10^-3 10^0 10^3 10^6])
legend(legstr)
hold on 
plot([-20 -20],[0 1],'--r')
hold on 
plot([20 20],[0 1],'--r')

%Model%%%%%%%%
subplot(2,3,6)
semilogx(IsomDom,pRodsModel,'k')
ylabel('rods/(rods+cones)')
xlabel('Isomerization Rate')
ylim([0 1])

set(gca,'XTick',IsomtickDom)
xlim([IsomtickDom(1) IsomtickDom(end)])

%CRT range
hold on
CRTlower = 0.8*10^4;
CRTupper = 0.8*10^4;
hold on
plot([CRTlower CRTlower],[0 1],'--b')
hold on 
plot([CRTupper CRTupper],[0 1],'--b')

%Solar range
solarlower = 10^1;
solarupper = 5*10^6;
hold on
plot([solarlower solarlower],[0 1],'--r')
hold on 
plot([solarupper solarupper],[0 1],'--r')

%Our display range
Projlower = 10000;
Projupper = 400000;
hold on
plot([Projlower Projlower],[0 1],'--k')
hold on 
plot([Projupper Projupper],[0 1],'--k')



%%



