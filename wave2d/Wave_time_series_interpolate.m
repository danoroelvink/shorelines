%% Convert form hourly to Daily wave data
clear all
close all
clc

WVC=load('Point_021821_ext.txt');

%col.1 (yyyymmddHHSS), col.2 (Hs), col.3 (Dir), col.4 (Tp)

i=1; %index of the original array
in=1; %index of the new array
% Hs=zeros(:length(WVC(:,1)));
n=336; %no. of data per year
Hs(:,1:ceil(length(WVC(:,1))/n))=0;
Ht(:,1:ceil(length(WVC(:,1))/n))=0;
Tr(:,1:ceil(length(WVC(:,1))/n))=0;
Hsin(:,1:ceil(length(WVC(:,1))/n))=0;
Hcos(:,1:ceil(length(WVC(:,1))/n))=0;
time(:,1:ceil(length(WVC(:,1))/n))=0;
x=round(length(WVC(:,1))/n)*n;

while i<x
    
    
    for ic=1:n
        
        time(in)=str2num(datestr(datenum(num2str(WVC(i,1)),'yyyymmddHHMM'),'yyyymmdd'));
        Hs(in)=Hs(in)+WVC(i,2)^2.5;
        Ht(in)=Ht(in)+WVC(i,2)^2;
        Tr(in)=Tr(in)+(WVC(i,2)^2)*WVC(i,3);
        Hsin(in)=Hsin(in)+(WVC(i,2)^2.5)*sind(WVC(i,4));
        Hcos(in)=Hcos(in)+(WVC(i,2)^2.5)*cosd(WVC(i,4));
        
        i=i+1;
    end
    
    
    in=in+1;
end

Hs=(Hs./n).^(1/2.5);
Tr=Tr./Ht;
Dirr=atan2d(Hsin,Hcos);

WVCn(:,1)=time;
WVCn(:,2)=Hs;
WVCn(:,3)=Tr;
WVCn(:,4)=Dirr;


% save('Damietta_daily_84to2010','WVCn');

%% Filtering out of range data


in=1; %index of the original array
inn=1; %index of the new array
% Hs=zeros(:length(WVC(:,1)));

%% if the data at this level is 0-360 use the lines below first 
c=length(find(WVCn(:,4) >=90 | WVCn(:,4) <= -90)); %no. of out of range values

xn=length(WVCn(:,1));
xnn=xn-c;   % length of new array

 WVCnn(1:xnn,1:4)=0;
while in<xn
    
    
    if WVCn(in,4) >=90 || WVCn(in,4)<= -90
        
        WVCnn(inn,1)=WVCn(in,1);
        WVCnn(inn,2)=0;
        WVCnn(inn,3)=0;
        WVCnn(inn,4)=0;
        in=in+1;
        inn=inn+1;
    else
        WVCnn(inn,1)=WVCn(in,1);
        WVCnn(inn,2)=WVCn(in,2);
        WVCnn(inn,3)=WVCn(in,3);
        WVCnn(inn,4)=WVCn(in,4);
        
        inn=inn+1;
        in=in+1;
    end
    
    
    
   
end

save('doubleWeekly_wave_table_delete','WVCnn','-ascii');

% figure(1)
% subplot(311)
% plot(WVC(:,2))
% subplot(312)
% plot(WVCn(:,2))
% subplot(313)
% plot(WVCnn(:,2))
% 
% figure(2)
% subplot(311)
% plot(WVC(:,3))
% subplot(312)
% plot(WVCn(:,3))
% subplot(313)
% plot(WVCnn(:,3))

figure(3)
subplot(311)
plot(WVC(:,4))
subplot(312)
plot(WVCn(:,4))
subplot(313)
plot(WVCnn(:,4))




%% additional in case of 360 data
% while in<xn
%     
%     
%     if WVCn(in,4) > 180
%        
%         WVCn(in,4)=WVCn(in,4)-360;
%         in=in+1;
%    
%     else
%       
%         WVCn(inn,4)=WVCn(in,4);
%         
%      
%         in=in+1;
%     end
%     
% 
%    
% end

