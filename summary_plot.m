clear

data = csv2cell("data.txt",";");

m1 = [data{:,3}];
m2 = [data{:,4}];
dist = [data{:,5}];
SNR = [data{:,6}];
RA = [data{:,7}];
dec = [data{:,8}];
iota = [data{:,9}];
area90 = [data{:,12}];
area50 = [data{:,13}];


threshold = 9
idx=find(SNR > threshold);

## luminosity distance (Mpc) redshift mass1 (Msun) mass2 (Msun) RA dec
load ../mergers.mat

clf

subplot(1,3,1)
hist(m1(idx))
xlabel("m_1 (Msun)");
title(sprintf("threshold SNR > %d -- %d binaries",threshold,length(idx)));
axis("square")

subplot(1,3,2)
[h1,n1]=hist([mergers.dist],'r;full population;');
h = plot(n1,h1,"r","linewidth",2); 
hold on
[h2,n2]=hist(dist(idx));
plot(n2,h2,"b","linewidth",2);
hold off
grid("on")
set(gca, 'yscale', 'log')'
xlabel("luminosity distance (Mpc)");
ylabel("count");
axis("square")
axis([0 500])
title("full pop (red) vs selected (blue)");

subplot(1,3,3)
hist(180/pi*iota(idx),'b')
hold off
grid("on")
xlabel("binary inclination (deg)");
ylabel("count");
axis("square")
axis([0 180])

print("-dpng","-FHelvetica:11","summary_plot.png");