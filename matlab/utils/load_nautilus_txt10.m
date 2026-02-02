function S = load_nautilus_txt10(filename)
%LOAD_NAUTILUS_TXT10 Read UST Nautilus IMU logs with 10 fields per row.
%
% Format (10 numeric fields per line, separated by whitespace or tabs):
%   seq  gx gy gz  ax ay az  mx my mz
%
% The magnetometer samples are returned in S.Y as a 3xN matrix.
% Gyro/accel fields are also parsed and returned for convenience, but the
% calibration solvers in this package use only S.Y.
%
% Output struct fields:
%   S.seq    Nx1 sequence numbers
%   S.gyro   3xN gyro samples
%   S.accel  3xN accel samples
%   S.mag    3xN magnetometer samples (same as S.Y)
%   S.Y      3xN magnetometer samples
%   S.N      number of samples
%   S.filename

fid = fopen(filename, 'r');
if fid < 0
    error('Could not open file: %s', filename);
end
c = onCleanup(@() fclose(fid));

% Robust whitespace parser (tabs or spaces). Also tolerate comma separators.
fmt = repmat('%f', 1, 10);
C = textscan(fid, fmt, 'Delimiter', {'\t', ' ', ','}, ...
    'MultipleDelimsAsOne', true, 'CollectOutput', true, 'CommentStyle', '#');

if isempty(C) || isempty(C{1})
    error('No numeric data found in %s', filename);
end

X = C{1};
if size(X,2) < 10
    error('Expected 10 columns but found %d in %s', size(X,2), filename);
end

seq = X(:,1);
gyro = X(:,2:4);
accel = X(:,5:7);
mag = X(:,8:10);

S = struct();
S.seq = seq;
S.gyro = gyro.';
S.accel = accel.';
S.mag = mag.';
S.Y = S.mag;
S.N = size(mag,1);
S.filename = filename;
end
