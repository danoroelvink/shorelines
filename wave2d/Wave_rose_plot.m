
x=load('wave_data_part1.txt');

wind_rose(x(:,4),x(:,2),'dtype','meteo','nAngles',16,'legStr','H[m]','titStr','Wave Height');
% print((strcat('Wave height')),'-dpng', '-r300')
%wind_rose(x(:,4),x(:,2),'dtype','meteo','nAngles',16,'legStr','T[sec]','titStr','peak period ');