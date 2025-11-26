% This script is used to emulate a jamming signal for a wireless communication competition.
% It continuously transmits complex white noise using the PlutoSDR device.

clc; clear; close all;

% PlutoSDR Connection
tx = sdrtx('Pluto');

% Transmission Parameters
tx.CenterFrequency = 2.4e9;
tx.BasebandSampleRate = 1e6;
tx.Gain = -10;
tx.ChannelMapping = 1;

% Generate Complex White Noise
fs = tx.BasebandSampleRate;
numSamples = fs * 1;
noise = (randn(numSamples,1) + 1i*randn(numSamples,1)) * 0.3;

% Continuous Transmission
transmitRepeat(tx, noise);

disp('Continuous white noise transmission started. Type "release(tx)" to stop.');
