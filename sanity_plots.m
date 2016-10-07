clear

data = csv2cell("data.txt",";");

m1 = [data{:,3}];
m2 = [data{:,4}];
dist = [data{:,5}];
SNR = [data{:,6}];
RA = [data{:,7}];
dec = [data{:,8}];
iota = [data{:,9}];

threshold = 10;
idx=find(SNR > threshold);

subplot(3,2,1)
hist(m1(idx))
xlabel("m1 Msun");

subplot(3,2,2)
hist(m2(idx))
xlabel("m2 Msun");

subplot(3,2,3)
hist(dist(idx))
xlabel("distance Mpc");

subplot(3,2,4)
hist(SNR(idx))
xlabel("SNR");
title(sprintf("SNR > %d -- %d binaries",threshold,length(idx)));

subplot(3,2,5)
hist(iota(idx))
xlabel("inclination deg");

subplot(3,2,6)
plot(RA(idx),dec(idx),".");
xlabel("RA");
ylabel("dec");

print("-dpng","-FHelvetica:12","sanity_plots.png");