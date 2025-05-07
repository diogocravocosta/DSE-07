clear

h_orbit = 600e3;
h_GS = 0;

dist = slantRangeCircularOrbit(5, h_orbit, h_GS)

Su = 2110e6;
Sd = 2290e6;
Xu = 7190e6;
Xd = 8500e6;

Kd = 27000e6;

frequ = [Su Xu];
freqd = [Sd Xd Kd];

fsplu = fspl(dist, freq2wavelen(frequ));
fspld = fspl(dist, freq2wavelen(freqd));

combinedu = 20 + fsplu - 228.6 + 50
combinedd = 20 + fspld - 228.6 + 70