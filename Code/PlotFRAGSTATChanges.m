%Plot vegetation change

%Once have values for all years, change this so all numbers in a matrix.

%load('VegAreaSlopeCctd.mat'); %Load areas and accuracy measures
load('C:/Users/boisrame/Google Drive/ThompsonFire/VegChangeAnals/FragstatsLandOut.mat'); %Load landscape level indices from Fragstats
load('C:/Users/boisrame/Google Drive/ThompsonFire/VegChangeAnals/FragstatClassOut.mat'); %Load class level metrics from Fragstats
load('VegPatches05m.mat'); %Load individual patch sizes from Fragstats
%load('FireList.mat'); %Load list of all fires
%load('FragstatClassOut2.mat'); %more class level metrics

%Year	 VegNum 	 PLAND 	 NP 	 LPI 	 FRAC_MN 	 FRAC_AM 	 FRAC_MD 	 FRAC_SD 	 FRAC_CV 	 PROX_MN 	 PROX_AM 	 PROX_MD 	 AI
YearC=1; VegNumC=2; LPIC=5; FracAMC=7; FRACsdC=9; ProxAMC=12; AggIndC=14;
YearL=1; SHEIL=2; SIEIL=3; AIL=4;
AreaMnC=5; AreaMdC=6; AreaRaC=7; AreaStdC=8;



Veg={'Shrub','Sparse Mdw.','Conifer','Dense Mdw.'};

color1=[.2 .65 .1];
color2=[.5 .3 .1];
color3=[.95 .95 .6];
color4=[.95 .65 0];

YearList=FragstatLand05m(:,YearL);

set(0,'defaultlinelinewidth',2)
set(0,'defaulttextfontsize',16);
set(0,'defaultaxesfontsize',14);

figure()
plot(FragstatLand05m(:,YearL),FragstatLand05m(:,SHEIL),'k--o');
hold on
plot(FragstatLand05m(:,YearL),FragstatLand05m(:,SIEIL),'k:s');
%plot(FragstatLand05m(:,YearL),FragstatLand05m(:,AIL),'k--s');
hold off
xlabel('Year');
ylabel('Index Value');
title('Landscape Evenness Indices','FontSize',16);
legend('SHEI','SIEI','Location','NorthWest');

figure()
plot(FragstatLand05m(:,YearL),FragstatLand05m(:,AIL),'k-.v');
hold off
xlabel('Year');
ylabel('Aggregation Index');
title('Landscape Aggregation Index');

figure()
subplot(2,1,2)
%plot(FragstatClass05m(1:5,YearL),FragstatClass05m(FragstatClass05m(:,VegNumC)==1,LPIC),'k-o');
hold on
plot(FragstatClass05m(1:5,YearL),FragstatClass05m(FragstatClass05m(:,VegNumC)==2,LPIC),'k--s');
plot(FragstatClass05m(1:5,YearL),FragstatClass05m(FragstatClass05m(:,VegNumC)==3,LPIC),'k:*');
plot(FragstatClass05m(1:5,YearL),FragstatClass05m(FragstatClass05m(:,VegNumC)==4,LPIC),'k-.v');
hold off
xlim([1965 2015]);
xlabel('Year');
ylabel('%');
legend(Veg{2:4},'Location','NorthWest');

subplot(2,1,1)
plot(FragstatClass05m(1:5,YearC),FragstatClass05m(FragstatClass05m(:,VegNumC)==1,LPIC),'k-o');
xlim([1965 2015]);xlabel('Year');
ylabel('%');
legend(Veg{1},'Location','SouthWest');
title('Largest Patch Percent Area');



%Mean patch area
figure()
subplot(2,1,2)
hold on
plot(FragstatClass05b(1:5,YearL),FragstatClass05b(FragstatClass05b(:,VegNumC)==2,AreaMnC),'k--s');
plot(FragstatClass05b(1:5,YearL),FragstatClass05b(FragstatClass05b(:,VegNumC)==3,AreaMnC),'k:*');
plot(FragstatClass05b(1:5,YearL),FragstatClass05b(FragstatClass05b(:,VegNumC)==4,AreaMnC),'k-.v');
hold off
xlim([1965 2015]);
xlabel('Year');
ylabel('ha');
%legend(Veg{2:4},'Location','NorthWest');

subplot(2,1,1)
plot(FragstatClass05b(1:5,YearC),FragstatClass05b(FragstatClass05b(:,VegNumC)==1,AreaMnC),'k-o');
xlim([1965 2015]);xlabel('Year');
ylabel('ha');
%legend(Veg{1},'Location','SouthWest');
title('Mean Patch Area');



figure()
subplot(2,1,2)
hold on
plot(FragstatClass05b(1:5,YearL),FragstatClass05b(FragstatClass05b(:,VegNumC)==2,AreaStdC),'k--s');
plot(FragstatClass05b(1:5,YearL),FragstatClass05b(FragstatClass05b(:,VegNumC)==3,AreaStdC),'k:*');
plot(FragstatClass05b(1:5,YearL),FragstatClass05b(FragstatClass05b(:,VegNumC)==4,AreaStdC),'k-.v');
hold off
xlim([1965 2015]);
xlabel('Year');
ylabel('ha');
%legend(Veg{2:4},'Location','NorthWest');

subplot(2,1,1)
plot(FragstatClass05b(1:5,YearC),FragstatClass05b(FragstatClass05b(:,VegNumC)==1,AreaStdC),'k-o');
xlim([1965 2015]);xlabel('Year');
ylabel('ha');
%legend(Veg{1},'Location','SouthWest');
title('Std. Dev. Patch Area');


%plot fire areas
FireYrs=unique(FireMat(:,1));
FireHa=0*FireYrs;
for i=1:length(FireYrs)
    FireHa(i)=sum(FireMat(FireMat(:,1)==FireYrs(i),2));
end
figure()
bar(FireYrs,FireHa,'FaceColor',[.2 .2 .2]);
hold on;
bar(YearList,3000*ones(size(YearList)),'FaceColor',[.7 .7 .7],'BarWidth',.001);
%plot([YearList(1) YearList(1)],[0 3000],'k--','LineWidth',1)
xlim([1965 2015]);
text(1966.5,2600,'1969'); text(1984.5,2600,'1987'); text(1994.5,2600,'1997'); text(2002.5,2600,'2005'); text(2009.5,2600,'2012');
xlabel('Year');
ylabel('Area Burned (Ha)');

TotPre87=sum(FireHa(FireYrs<1987))
Tot_87_97=sum(FireHa(FireYrs<1997&FireYrs>=1987))
Tot_97_05=sum(FireHa(FireYrs<2005&FireYrs>=1997))
Tot_05_12=sum(FireHa(FireYrs<2012&FireYrs>=2005))

figure()
subplot(2,1,2)
bar(FireYrs,FireHa);
hold on
plot(FragstatClass05m(1:5,YearL),FragstatClass05m(FragstatClass05m(:,VegNumC)==2,LPIC),'k--s');
plot(FragstatClass05m(1:5,YearL),FragstatClass05m(FragstatClass05m(:,VegNumC)==3,LPIC),'k:*');
plot(FragstatClass05m(1:5,YearL),FragstatClass05m(FragstatClass05m(:,VegNumC)==4,LPIC),'k-.v');
hold off
xlabel('Year');
ylabel('%');
legend(Veg{2:4},'Location','NorthWest');

subplot(2,1,1)
plot(FragstatClass05m(1:5,YearC),FragstatClass05m(FragstatClass05m(:,VegNumC)==1,LPIC),'k-o');
xlabel('Year');
ylabel('%');
legend(Veg{1},'Location','SouthWest');
title('Largest Patch Percent Area');

VegPatches(:,4)=VegPatches(:,3);
VegPatches(VegPatches(:,3)>5,4)=5; %Set max to 20ha so can see distribution
VegPatches(VegPatches(:,3)<.01,4)=NaN;
clnum=1; %Which Class
coln=4; %Use 4 if want to truncate at a certain size
ymx=300;
%xmx=2000;
nb=100;

h69 = histogram(VegPatches((VegPatches(:,1)==1969)&(VegPatches(:,2)==clnum),coln),nb);
h12 = histogram(VegPatches((VegPatches(:,1)==2012)&(VegPatches(:,2)==clnum),coln),nb);

figure()
semilogx(h69.BinEdges(2:end),h69.Values,'g')
hold on;
semilogx(h12.BinEdges(2:end),h12.Values,'k')
hold off;

h69l = histogram(log10(VegPatches((VegPatches(:,1)==1969)&(VegPatches(:,2)==clnum),coln)),nb);
h12l = histogram(log10(VegPatches((VegPatches(:,1)==2012)&(VegPatches(:,2)==clnum),coln)),nb);
figure()
plot(h69l.BinEdges(2:end),h69.Values,'g')
hold on;
plot(h12l.BinEdges(2:end),h12.Values,'k')
hold off;


figure()
subplot(3,2,1)
ylim([0 ymx])
title('Sparse Patch Sizes - 1969');
xlabel('Area (ha)');
ylabel('Number');
subplot(3,2,2)
h87=histogram(VegPatches((VegPatches(:,1)==1987)&(VegPatches(:,2)==clnum),coln),nb);
ylim([0 ymx]); title('1987');
xlabel('Area (ha)');
ylabel('Number');
subplot(3,2,3)
histogram(VegPatches((VegPatches(:,1)==1997)&(VegPatches(:,2)==clnum),coln),nb);
ylim([0 ymx]); title('1997');
xlabel('Area (ha)');
ylabel('Number');
subplot(3,2,4)
h05=histogram(VegPatches((VegPatches(:,1)==2005)&(VegPatches(:,2)==clnum),coln),nb);
ylim([0 ymx]); title('2005');
xlabel('Area (ha)');
ylabel('Number');
subplot(3,2,5)
ylim([0 ymx]); title('2012');
xlabel('Area (ha)');
ylabel('Number');




%A power law relationship
%[slope, intercept] = logfit(x,y,'loglog');
%           yApprox = (10^intercept)*x.^(slope);
x=h69.BinEdges(2:end);
y=h69.Values;
FL = fit(log10(x(y>0)'),log10(y(y>0)'),'poly1');
slp=FL.p1; int=FL.p2;
yfit=(10^int)*x.^slp;
corr(yfit',y')

x=h12.BinEdges(2:end);
y=h12.Values;
FL = fit(log10(x(y>0)'),log10(y(y>0)'),'poly1');
slp=FL.p1; int=FL.p2;
yfit=(10^int)*x.^slp;
corr(yfit',y')

%All class types:
h69a = histogram(VegPatches(VegPatches(:,1)==1969,3),100);
 Max69 = max(VegPatches(VegPatches(:,1)==1969,3))
 Max12 = max(VegPatches(VegPatches(:,1)==2012,3))
 P25_69 = prctile(VegPatches(VegPatches(:,1)==1969,3),25)
 P25_12 = prctile(VegPatches(VegPatches(:,1)==2012,3),25)
 P50_69 = prctile(VegPatches(VegPatches(:,1)==1969,3),50)
 P50_12 = prctile(VegPatches(VegPatches(:,1)==2012,3),50)

x=h69a.BinEdges(2:end);
y=h69a.Values;
FL = fit(log10(x(y>0)'),log10(y(y>0)'),'poly1');
slp=FL.p1; int=FL.p2;
yfit=(10^int)*x.^slp;
[r,p]=corr(yfit',y')

h12a = histogram(VegPatches(VegPatches(:,1)==2012,3),100);
x=h12a.BinEdges(2:end);
y=h12a.Values;
FL = fit(log10(x(y>0)'),log10(y(y>0)'),'poly1');
slp=FL.p1; int=FL.p2;
yfit=(10^int)*x.^slp;
[r,p]=corr(yfit',y')


FracMat=FragstatClass05m(:,[YearC VegNumC FracAMC FRACsdC]);
FracMat=FracMat(FracMat(:,1)~=1997,:);
figure()
plot(FracMat(1:4,1),FracMat(FracMat(:,2)==1,3),'k-o');
hold on
plot(FracMat(1:4,1),FracMat(FracMat(:,2)==2,3),'k--s');
plot(FracMat(1:4,1),FracMat(FracMat(:,2)==3,3),'k:*');
plot(FracMat(1:4,1),FracMat(FracMat(:,2)==4,3),'k-.v');
%errorbar(FracMat(:,1),FracMat(:,3),FracMat(:,4));
hold off
xlabel('Year')
ylabel('Fractal Dimension');
title('Area-Weighted Fractal Dimension');
legend(Veg,'Location','West');


figure()
plot(FragstatClass05m(1:5,1),FragstatClass05m(FragstatClass05m(:,VegNumC)==1,AggIndC),'k-o');
hold on
plot(FragstatClass05m(1:5,1),FragstatClass05m(FragstatClass05m(:,VegNumC)==2,AggIndC),'k--s');
plot(FragstatClass05m(1:5,1),FragstatClass05m(FragstatClass05m(:,VegNumC)==3,AggIndC),'k:*');
plot(FragstatClass05m(1:5,1),FragstatClass05m(FragstatClass05m(:,VegNumC)==4,AggIndC),'k-.v');
hold off
xlabel('Year')
ylabel('Aggregation Index');
title('Class Aggregation Index');
legend(Veg,'Location','NorthEast');

figure()
subplot(2,1,1);
plot(FragstatClass05m(1:5,1),FragstatClass05m(FragstatClass05m(:,VegNumC)==1,ProxAMC),'k-o');
title('Class Proximity Index');
ylabel('Proximity Index');
legend(Veg{1},'Location','NorthWest')
subplot(2,1,2);
hold on
plot(FragstatClass05m(1:5,1),FragstatClass05m(FragstatClass05m(:,VegNumC)==2,ProxAMC),'k--s');
plot(FragstatClass05m(1:5,1),FragstatClass05m(FragstatClass05m(:,VegNumC)==3,ProxAMC),'k:*');
plot(FragstatClass05m(1:5,1),FragstatClass05m(FragstatClass05m(:,VegNumC)==4,ProxAMC),'k-.v');
hold off
xlabel('Year')
ylabel('Proximity Index');
legend(Veg{2:4},'Location','NorthWest');

ProxNormVec=reshape(FragstatClass05m(:,ProxAMC),5,4);
ProxNormVec=ProxNormVec./repmat(max(ProxNormVec),5,1);
figure()
hold on;
plot(FragstatClass05m(1:5,1),ProxNormVec(:,1),'k-o');
plot(FragstatClass05m(1:5,1),ProxNormVec(:,2),'k--s');
plot(FragstatClass05m(1:5,1),ProxNormVec(:,3),'k:*');
plot(FragstatClass05m(1:5,1),ProxNormVec(:,4),'k-.v');
hold off;
title('Normalized Proximity Index');
xlabel('Year');
legend(Veg,'Location','SouthEast');

%Vegetation Total Area from Maps
VegArea_km2 = VegArea_m2*4/1000000; %Multiply by 4 to deal w/ error in area calculated before (used 25m2 grid when actually was 100m2)
VegA_e = (1-VegArea_Rel).*VegArea_km2;


figure()
subplot(3,1,1)
grid on;
title('Total Area');
plot(YearList,VegArea_km2(:,1),'k-o'); hold on;
errorbar(YearList,VegArea_km2(:,1),VegA_e(:,1),'k.');
title('Total Area');
ylabel('km^2')
legend(Veg{1},'Location','SouthWest')
subplot(3,1,2)
plot(YearList,VegArea_km2(:,2),'k--s',YearList,VegArea_km2(:,3),'k:*');
hold on;
errorbar(YearList,VegArea_km2(:,2),VegA_e(:,2),'k.');
errorbar(YearList,VegArea_km2(:,3),VegA_e(:,3),'k.');
ylabel('km^2')
legend(Veg{2:3},'Location','NorthWest')
subplot(3,1,3)
plot(YearList,VegArea_km2(:,4),'k-.v');
hold on;
errorbar(YearList,VegArea_km2(:,4),VegA_e(:,4),'k.','MarkerSize',1);
legend(Veg{4},'Location','NorthWest')
xlabel('Year');
ylabel('km^2')


figure()
b=bar(VegArea_km2');
b5c=(2025-2012)/(2025-1969);
set(b(5),'FaceColor',[b5c b5c b5c]);
b4c=(2025-2005)/(2025-1969);
set(b(4),'FaceColor',[b4c b4c b4c]);
b3c=(2025-1997)/(2025-1969);
set(b(3),'FaceColor',[b3c b3c b3c]);
b2c=(2025-1987)/(2025-1969);
set(b(2),'FaceColor',[b2c b2c b2c]);
b1c=(2025-1969)/(2025-1969);
set(b(1),'FaceColor',[b1c b1c b1c]);
set(b,'EdgeColor',[0 0 0])
hold on;
errorbar(repmat(1:4,5,1)+repmat([-.3 -.15 0 .15 .3]',1,4),VegArea_km2,VegA_e,'k.','MarkerSize',1);
l1=legend('1969','1987','1997','2005','2012','Location','NorthEast');
h=gca;
h.XTickLabel=Veg;
h.FontSize=13;
set(l1,'FontSize',14);
xlabel('Vegetation Class');
ylabel('km^2');
title('Total Area');

%Fix up area of 1997 map being a little lower in total area
VegArea_km2b=VegArea_km2;
VegArea_km2b(3,1:3)=VegArea_km2(3,1:3)+repmat(0.8,1,3);
figure()
b=bar([1969 1987 1997 2005 2012], VegArea_km2b,'stacked');
%Black and White
%set(b(1),'FaceColor',[.95 .95 .95])
%set(b(2),'FaceColor',[.7 .7 .7])
%set(b(3),'FaceColor',[.4 .4 .4])
%set(b(4),'FaceColor',[0 0 0])
%Color
set(b(1),'FaceColor',[.05 .65 .05])
set(b(2),'FaceColor',[.7 .4 .2])
set(b(3),'FaceColor',[.8 .8 0])
set(b(4),'FaceColor',[.6 .4 0])
xlabel('Year')
ylabel('Area (km^2)')
ylim=[0,110];
legend('Conifer','Shrub','Sparse Meadow','Dense Meadow','Location','EastOutside');

f1=figure();
ax1 = axes('Parent',f1,'XTickLabel',{'1969','1987','1997','2005','2012'},'XTick',[1 2 3 4 5]);
box(ax1,'on');
hold(ax1,'on');
b=bar(VegArea_km2b,'stacked');
%set(b(1),'FaceColor',[.95 .95 .95])
%set(b(2),'FaceColor',[.7 .7 .7])
%set(b(3),'FaceColor',[.4 .4 .4])
%set(b(4),'FaceColor',[0 0 0])
set(b(1),'FaceColor',[.2 .65 .1])
set(b(2),'FaceColor',[.5 .3 .1])
set(b(3),'FaceColor',[.95 .95 .6])
set(b(4),'FaceColor',[.95 .65 0])
xlabel('Year')
ylabel('Area (km^2)')
ylim=[0,110];
legend('Conifer','Shrub','Sparse Meadow','Dense Meadow','Location','EastOutside');


%Changes
YearlyChange=100*(VegArea_m2(2:5,:)-VegArea_m2(1:4,:))./VegArea_m2(1:4,:)
YearlyChangeRate=YearlyChange./(repmat([1987-1973;1997-1987;2005-1997;2012-2005],1,4));
TotChange=100*(VegArea_m2(5,:)-VegArea_m2(1,:))./(VegArea_m2(1,:))
PerYear1=mean(YearlyChangeRate)
PerYear2=median(YearlyChangeRate)
PerYear3=TotChange/(2012-1973)

BurnGrouped=[TotPre87 Tot_87_97 Tot_97_05 Tot_05_12]';
PerBurnArea=(VegArea_m2(2:5,:)-VegArea_m2(1:4,:))./repmat(BurnGrouped,1,4)



%load('VegChangeVals2.mat');
%load('LandIndices.mat');

%Year69=Year69/1000000; %Convert area to km^2
%Year12=Year12/1000000;

% Veg={'Conifer','Shrub','Sp. Grass','Wet Meadow','Aspen'};
% 
% Year69=Year69(1:4);
% Year12=Year12(1:4);
% Veg=Veg(1:4);
% Rel69=Rel69(1:4);
% Rel97=Rel97(1:4);
% Rel12=Rel12(1:4);
% 
% 
% Diff_69_12=100*(Year12-Year69)./Year69;
% 
% Err69=(1-Rel69);
% Err97=(1-Rel97);
% Err12=(1-Rel12);
% 
% DU=((Year12.*(1+Err12))-Year69.*(1-Err69))./Year69.*(1-Err69)
% DL=((Year12.*(1-Err12))-Year69.*(1+Err69))./Year69.*(1+Err69)
% DM1=((Year12.*(1+Err12))-Year69.*(1+Err69))./Year69.*(1+Err69)
% DM2=((Year12.*(1-Err12))-Year69.*(1-Err69))./Year69.*(1-Err69)
% De1=max([DU,DL,DM1,DM2]')*100-Diff_69_12'
% De2=Diff_69_12'-min([DU,DL,DM1,DM2]')*100
% 
% Year97=Year69; %TEMPORARY!!!!!!!!!!!!!!
% figure()
% b=bar([Year69,Year97,Year12],1)
% b3c=(2020-2012)/(2012-1969);
% set(b(3),'FaceColor',[b3c b3c b3c]);
% b2c=(2020-1997)/(2020-1969);
% set(b(2),'FaceColor',[b2c b2c b2c]);
% set(b(1),'FaceColor',[1 1 1]);
% hold on
% errorbar([(1:4)-.2;1:4;(1:4)+.2]',[Year69,Year97,Year12],[Err69.*Year69, Err97.*Year97, Err12.*Year12],'k.')
% h=gca;
% h.XTickLabel=Veg;
% h.FontSize=13;
% %set(h.XTickLabel,'FontSize',16)
% l1=legend('1969','1997','2012');
% set(l1,'FontSize',14);
% %xlabel('Vegetation Type','FontSize',14)
% ylabel('Total Area (km^2)','FontSize',14)
% 
% 
% Res=LI(:,1);
% Years=LI(1:5,2);
% figure()
% plot(Years,LI(Res==5,3))
% hold on
% plot(Years,LI(Res==10,3));
% hold off;
% ylabel('SHEI')
% legend('5m','10m');
% title('Shannons Evenness Index');
% 
% Res=LI(:,1);
% Years=LI(1:5,2);
% figure()
% plot(Years,LI(Res==5,4))
% hold on
% plot(Years,LI(Res==10,4));
% hold off;
% ylabel('SIEI');
% legend('5m','10m');
% 
% Res=LI(:,1);
% Years=LI(1:5,2);
% figure()
% plot(Years,LI(Res==5,5))
% hold on
% plot(Years,LI(Res==10,5));
% hold off;
% ylabel('AI');
% legend('5m','10m');

%CHANGES 1969-2012 ONLY

% figure()
% bar(1:4,Diff_69_12,'w')
% hold on
% errorbar(1:4,Diff_69_12,De1,De2,'k.')
% %xlabel('Vegetation Type')
% ylabel('% Change in Area','FontSize',14)
% h=gca;
% h.XTickLabel=Veg;
% h.FontSize=13;

% figure()
% b=bar([Year69,Year12],1)
% set(b(2),...
%     'FaceColor',[0.501960784313725 0.501960784313725 0.501960784313725]);
% set(b(1),'FaceColor',[1 1 1]);
% hold on
% errorbar([(1:4)-.15;(1:4)+.15]',[Year69,Year12],[Err69.*Year69, Err12.*Year12],'k.')
% h=gca;
% h.XTickLabel=Veg;
% h.FontSize=13;
% %set(h.XTickLabel,'FontSize',16)
% l1=legend('1969','2012');
% set(l1,'FontSize',14);
% %xlabel('Vegetation Type','FontSize',14)
% ylabel('Total Area (km^2)','FontSize',14)
% 
% figure()
% b=bar([1969,2012],[Year69,Year12]',.8,'stacked');
% l1=legend(Veg);
% ylim([0,132.3])
% set(l1,'FontSize',14);
% ylabel('Total Area (km^2)','FontSize',14);
% 
% figure()
% b=bar([1969,2012],[Year69/sum(Year69),Year12/sum(Year12)]',.8,'stacked');
% ylim([0,1])
% l1=legend(Veg);
% set(l1,'FontSize',14);
% ylabel('Proportional Area');
