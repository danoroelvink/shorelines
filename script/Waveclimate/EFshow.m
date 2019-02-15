clear all; close all;

Input = 'wavedata_IJMDMNTSPS_200501010000_201412312359.txt';
%Input = 'wavedata_SCHIERMNOND_200101010000_201012312359.txt';
Data2  = load(Input);
maxi  = size(Data2,1);
rho   = 1025;           %%water density
g     = 9.81;           %%gravity acceleration
MF    = 100

dirbin=50;              
hsbin=1;

hn=4;
mindeldir=0;            %% leave it for backward compatibility
maxdeldir=0;
datalowerlimit = 200;
dataupperlimit = 360;

jadefinedirbins = 0;    %% switch to define wave directions yourself (1) or not (0)
dirbinlims = [ 260 315 360];

Data = Data2(Data2(:,4)>datalowerlimit & Data2(:,4)<dataupperlimit,:);

table(:,1)=Data(:,2);                     % Hs
table(:,2)=Data(:,4);                     % Dir 
table(table(:,2)==0,2) = 360;
table(:,3)=Data(:,3);                     % Tp
table(:,4)=(rho*g*table(:,1).^2/8);       % energy
table(:,4)=1.56*table(:,4).*table(:,3)/2; % energy times Cg = flux

r=sortrows(table,2);                      % sort by directions

l1  = find(r(:,2)<mindeldir,1,'last');
l2  = find(r(:,2)>maxdeldir,1','first');

table2 = r(1:l1,:);
table3 = r(l2:size(r),:);

sizen               = (size(table2)+size(table3));
tablen(1:sizen,1:4) = 0;

tablen(1:size(table3),:) = table3;
tablen(size(table3)+1:size(table3)+size(table2),:) = table2;

clear r;
r    = tablen;
maxi = size(r,1);      %% no of wave points
SE   = sum(r(:,4));    %% total energy flux
j    = 0;
i    = 0;
% disp(['after first run: ',num2str(SE)]);
if ~jadefinedirbins
    for k = 1:dirbin
       AE = 0;
       while AE < (SE/dirbin)
          i = i+1;
          if i > maxi
             i = i-1;
             %break
          end
          AE = r(i,4)+AE;      %% sum until 1/dirbin of energy flux
       end
       D(k) = r(i,2);          %% direction limit
       l(k) = i;               %% last observation number
    end
else
   dirbin = length(dirbinlims);
   for k=1:dirbin
      i    = find(r(:,2)<dirbinlims(k),1,'last');
      D(k) = r(i,2);
      l(k) = i;
   end
end

%% Hs bins
% sort the data in the directional bins along wave heights
temp = tablen(1:l(1),:);  %%  all the observations in this bin
temp = sortrows(temp,1);  %%  sort on wave height
tablen(1:l(1),:) = temp;  %% rewrite in table
clear temp
for i=1:size(D')-1      
   temp = tablen(l(i)+1:l(i+1),:);
   temp = sortrows(temp,1);
   tablen(l(i)+1:l(i+1),:) = temp; 
   clear temp
end

clear r
r    = tablen;
maxi = size(r,1);
SE   = sum(r(:,4));
% disp(['after second run: ',num2str(SE)]);
% assign wave height to classes
if ~jadefinedirbins
    j=0;
    i=0;
    for k = 1:dirbin*hsbin
       AE = 0;
       while (AE<(SE/(dirbin*hsbin)) && i<l(floor(k/hsbin-0.00001)+1))
          i=i+1;
          if i>maxi
             i=i-1;
          end
          AE=r(i,4)+AE;
       end
       HH(k)=r(i,1);     %% wave height of last observation point in class
       lH(k)=i;          %% index last observation
    end
else
    clear SE;
    SE = zeros(dirbin,1);
    SE(1) = sum(r(1:l(1),4));
    for k=2:dirbin
       SE(k) = sum(r(l(k-1):l(k),4));    %% partial sums of energy flux 
    end
    i = 0;
    for k=1:dirbin
         ihs = 0;
         AE = 0;
         for ihs = 1:hsbin
            ll = (k-1)*hsbin + ihs;
            while AE < SE(k)*ihs/hsbin && i<l(floor(ll/hsbin-0.00001)+1)
               i = i+1;
               if i>maxi
                  i=i-1;
               end
               AE=r(i,4)+AE;
            end
            HH(ll) = r(i,1);
            lH(ll) = i;
         end
    end
end
%% Rep. values
binlim(2:dirbin*hsbin+1)=lH;
binlim(1)=1;
% Wavecase 1:Hs mean - 2:Dir mean - 3:Tp mean - 4:Sum Ef - 5:Hs rep based on
% Ef - 6:number of waves - 7:probability of occ.
% 
for i=1 : size(binlim')-1
    wcase(i,1)= mean(tablen(binlim(i):binlim(i+1),1));
    cosine(i) = mean(cos(deg2rad(tablen(binlim(i):binlim(i+1),2))));
    sine(i)   = mean(sin(deg2rad(tablen(binlim(i):binlim(i+1),2))));
    wcase(i,2)= mod(rad2deg(atan2(sine(i),cosine(i))),360);
    wcase(i,3)= mean(tablen(binlim(i):binlim(i+1),3));
    wcase(i,4)= sum(tablen(binlim(i):binlim(i+1),4));
    cg=1.56*wcase(i,3)/2;
    wcase(i,5)=((wcase(i,4)/(binlim(i+1)-binlim(i)))/cg/rho/g*8)^0.5;
    wcase(i,6)=(binlim(i+1)-binlim(i));
end

wcase(:,7)=wcase(:,6)/sum(wcase(:,6))*100*length(Data)/length(Data2);

% Plot
figure(1);
set(gcf,'color', 'w')
plot(tablen(:,2),tablen(:,1),'.');
set(gca,'fontsize',12)
grid on
hold on
HPL(1:dirbin,1:hsbin+1)=0;
for i=1 : dirbin
    HPL(i,2:hsbin+1)=HH(hsbin*(i-1)+1:hsbin*i);
end
DP(2:dirbin+1)=D;
DP(1)=min(tablen(:,2));

for i=1 : dirbin
    plot ([DP(i) DP(i)],[HPL(i,1) HPL(i,hsbin+1)],'r-')
    plot ([DP(i+1) DP(i+1)],[HPL(i,1) HPL(i,hsbin+1)],'r-')
    for j=1 : hsbin+1
        plot([DP(i) DP(i+1)],[HPL(i,j) HPL(i,j)],'r-')
    end
end
xlabel('wave directions [deg from N]','fontsize',12);
ylabel('wave heights [m]','fontsize',12);
title('');

scatter(wcase(:,2),wcase(:,5),'ro','filled')
%scatter(wcase(:,2),wcase(:,1),'g*')

wcon(:,1)=wcase (:,5);
wcon(:,2)=wcase (:,3);
wcon(:,3)=wcase (:,2);
wcon(:,4)=wcase (:,7)/100*365;

save ('wcon.mat','wcon') 
save ('wcon.txt','wcon','-ascii') 
