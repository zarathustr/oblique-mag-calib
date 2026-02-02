function data = load_nautilus_dataset(folder, pattern)
%LOAD_NAUTILUS_DATASET Load and concatenate UST Nautilus TXT logs.
%
% data = load_nautilus_dataset(folder)
% data = load_nautilus_dataset(folder, pattern)
%
% Inputs:
%   folder  path containing one or more .txt logs
%   pattern optional filename pattern, default 'ust-nautilus-mag-calib-*.txt'
%
% Output struct:
%   data.Y      3xN concatenated magnetometer samples
%   data.seq    1xN concatenated sequence numbers
%   data.files  cell array of file paths (sorted)
%   data.N      total number of samples
%
% Notes:
% - Each file is parsed with load_nautilus_txt10.
% - Sequence numbers are concatenated in arrival order; no attempt is made
%   to resample or align clocks across different sessions.

if nargin < 2 || isempty(pattern)
    pattern = 'ust-nautilus-mag-calib-*.txt';
end

if ~isfolder(folder)
    error('Folder does not exist: %s', folder);
end

L = dir(fullfile(folder, pattern));
if isempty(L)
    error('No files matched pattern %s in %s', pattern, folder);
end

% Sort by name for reproducibility
[~, idx] = sort({L.name});
L = L(idx);

files = cell(1, numel(L));
Yall = zeros(3, 0);
seqall = zeros(1, 0);

for i = 1:numel(L)
    f = fullfile(folder, L(i).name);
    files{i} = f;
    S = load_nautilus_txt10(f);
    Yall = [Yall, S.Y]; %#ok<AGROW>
    seqall = [seqall, S.seq(:).']; %#ok<AGROW>
end

data = struct();
data.Y = Yall;
data.seq = seqall;
data.files = files;
data.N = size(Yall,2);
end
