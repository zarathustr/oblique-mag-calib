function data = load_0711_dataset(folder)
%LOAD_0711_DATASET Load the attached 0711 magnetometer dataset.
%
% data = load_0711_dataset(folder)
%
% Inputs:
%   folder  directory containing 1.csv,2.csv,...
%
% Output struct:
%   data.Y        3xN concatenated magnetometer samples
%   data.N        total number of samples
%   data.files    cell array of file paths
%   data.N_per    samples per file

if nargin < 1 || isempty(folder)
    folder = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', '0711');
end

csvs = dir(fullfile(folder, '*.csv'));
if isempty(csvs)
    error('No CSV files found in folder: %s', folder);
end

% Sort by numeric prefix if possible
names = {csvs.name};
nums = nan(size(names));
for i = 1:numel(names)
    tok = regexp(names{i}, '^(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        nums(i) = str2double(tok{1});
    end
end
[~, order] = sortrows([isnan(nums(:)), nums(:)]);
csvs = csvs(order);

Yall = zeros(3, 0);
N_per = zeros(numel(csvs), 1);
files = cell(numel(csvs), 1);
for i = 1:numel(csvs)
    f = fullfile(folder, csvs(i).name);
    S = load_3dmg_mag_csv(f);
    Yall = [Yall, S.Y]; %#ok<AGROW>
    N_per(i) = S.N;
    files{i} = f;
end

data = struct();
data.Y = Yall;
data.N = size(Yall, 2);
data.files = files;
data.N_per = N_per;
data.folder = folder;
end
