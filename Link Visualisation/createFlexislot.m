function ant = createFlexislot(max_frequency)
%CREATEFLEXISLOT Summary of this function goes here
%   Detailed explanation goes here

%% Declare azimuth gain pattern
azangles = 0:5:360;
azangles_fine = 0:1:360;
azpattern_dB = [5 5 5 5 4.5 4 3.5 3 2.5 2 1.5 1 0 -0.75 -1.5 -2.25 -3 -4 -5 -6.25 -7.5 -9 -11 -13 -15 -18 -22 -23.5 -24.5 -24 -23.5 -22.5 -21 -20 -19.5 -19.25 -19];
azpattern_dB = [azpattern_dB flip(azpattern_dB(1:end-1))];
azpattern_dB_fine = interp1(azangles, azpattern_dB, azangles_fine, "spline");
%polarplot(deg2rad(azangles), azpattern_dB)
%rlim([-35 6])

%% Declare elevation gain pattern
elangles = -90:5:90;
elangles_fine = -90:1:90;
elpattern_dB = [4.5 4.5 5 5.5 5 4.5 6 5.5 5 6.5 6.5 5.5 6 7 7 7 6.5 5.5 4];
elpattern_dB = [flip(elpattern_dB) elpattern_dB(2:end)];
elpattern_dB_fine = interp1(elangles, elpattern_dB, elangles_fine, "spline");
%polarplot(deg2rad(elangles), elpattern_dB)
%rlim([-35 9])

%% Create 3D pattern from slices
thetaangles = 90 - elangles_fine;
[antpattern, patterntheta, patternphi] = patternFromSlices(elpattern_dB_fine, thetaangles, azpattern_dB_fine, azangles_fine, "Method", "CrossWeighted");
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

