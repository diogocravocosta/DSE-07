clear
clc

%% Define scenario
startTime = datetime(2032,1,1,12,00,0);      % 1 January 2032 12:00 UTC
stopTime = startTime + hours(1);             % 4 January 2032 12:00 UTC
sampleTime = 5;                              % seconds

scenario = satelliteScenario(startTime,stopTime,sampleTime);

%% Add vehicle
h2ermes = vehicleInOrbit(scenario);

%% Add ground station
gs = createGroundStation(scenario);

pointAt(gs.Gimbals, h2ermes);

%% Play scenario and visualise links

uplink = link(gs.Gimbals.Transmitters, h2ermes.Receivers);
downlink = link(h2ermes.Transmitters, gs.Gimbals.Receivers);

pattern(h2ermes.Transmitters, Size=100);
pattern(gs.Gimbals.Transmitters, Size=300000);

play(scenario)

%% Plot margin for uplink and downlink

[uplink_e, uplink_t] = ebno(uplink);
uplink_margin = uplink_e - h2ermes.Receivers.RequiredEbNo(1);
plot(uplink_t, uplink_margin);
xlabel("Time");
ylabel("Link margin (dB)");
title("Uplink margin");
grid on;

figure;
[downlink_e, downlink_t] = ebno(downlink);
downlink_margin = max(downlink_e, [], 1) - gs.Gimbals.Receivers.RequiredEbNo;
plot(downlink_t, downlink_margin);
xlabel("Time");
ylabel("Link margin (dB");
title("Downlink margin");
grid on;
