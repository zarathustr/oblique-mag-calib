function S = load_3dmg_mag_csv(filename)
%LOAD_3DMG_MAG_CSV Robustly read MicroStrain 3DM-GX5 style CSV logs.
%
% The 0711 dataset contains metadata blocks, then a header line starting with
% "GPS TFlags" and a sequence of comma-separated samples. This function
% searches for that header and extracts only the magnetometer channels
% (X Mag, Y Mag, Z Mag). It tolerates extra columns and occasional corrupted
% header tokens by selecting columns via substring match.
%
% Output:
%   S.Y        3xN magnetometer samples
%   S.N        number of samples
%   S.filename source filename

fid = fopen(filename, 'r');
if fid < 0
    error('Could not open file: %s', filename);
end
c = onCleanup(@() fclose(fid));

hdr = '';
while true
    tline = fgetl(fid);
    if ~ischar(tline)
        error('No data header line found in %s', filename);
    end
    if startsWith(tline, 'GPS TFlags')
        hdr = tline;
        break;
    end
end

cols = strsplit(hdr, ',');
ix = find(contains(cols, 'X Mag'), 1, 'first');
iy = find(contains(cols, 'Y Mag'), 1, 'first');
iz = find(contains(cols, 'Z Mag'), 1, 'first');
if isempty(ix) || isempty(iy) || isempty(iz)
    error('Mag columns not found in %s', filename);
end
need = max([ix, iy, iz]);

x = zeros(1, 0);
y = zeros(1, 0);
z = zeros(1, 0);

while true
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end
    tline = strtrim(tline);
    if isempty(tline)
        continue;
    end
    parts = strsplit(tline, ',');
    if numel(parts) < need
        continue;
    end
    vx = str2double(parts{ix});
    vy = str2double(parts{iy});
    vz = str2double(parts{iz});
    if any(isnan([vx, vy, vz]))
        continue;
    end
    x(end+1) = vx; %#ok<AGROW>
    y(end+1) = vy; %#ok<AGROW>
    z(end+1) = vz; %#ok<AGROW>
end

S = struct();
S.Y = [x; y; z];
S.N = numel(x);
S.filename = filename;
end
