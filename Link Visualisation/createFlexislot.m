function ant = createFlexislot(max_frequency)
%CREATEFLEXISLOT Summary of this function goes here
%   Detailed explanation goes here

%% Declare azimuth gain pattern
azangles = 0:10:360;
azpattern_dB = [5 5 4 4 3 2 0 -2 -4 -6 -8 -10 -14 -20 -24 -23 -21 -19 -18];
azpattern_dB = [azpattern_dB flip(azpattern_dB(1:end-1))];
%polarplot(deg2rad(azangles), azpattern_dB)
%rlim([-35 6])

%% Declare elevation gain pattern
elangles = -90:10:90;
%elpattern_dB = [5 5 5 6 5 7 7 8 7.5 5 -1 -6 -4 -8 -11 -14 -16 -17 -18];
elpattern_dB = [5 5 5 6 5 7 7 8 7.5 5];
elpattern_dB = [flip(elpattern_dB) elpattern_dB(2:end)];
%polarplot(deg2rad(elangles), elpattern_dB)
%rlim([-35 9])

%% Create 3D pattern from slices
thetaangles = 90 - elangles;
[antpattern, patterntheta, patternphi] = patternFromSlices(elpattern_dB, thetaangles, azpattern_dB, azangles, "Method", "CrossWeighted");
%patternFromSlices(elpattern_dB, thetaangles, azpattern_dB, azangles, "Method", "CrossWeighted")
%figure;
%patternFromSlices(elpattern_dB, thetaangles, azpattern_dB, azangles, "Method", "Summing")

%% Create custom antenna with given pattern
ant = phased.CustomAntennaElement("PatternCoordinateSystem", "phi-theta", ...
    "FrequencyVector", [0 max_frequency], ...
    "ThetaAngles", flip(patterntheta), ...
    "PhiAngles", patternphi, ...
    "MagnitudePattern", antpattern.', ...
    "PhasePattern", zeros(size(antpattern.')));
end

