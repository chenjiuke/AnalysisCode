function S = readJuliaDSC(fpath)
%READJULIADSC  Read a JULIA / Anton Paar DSC export (tab-separated, '#' header).
%
%   S = readJuliaDSC(fpath) returns a struct with fields
%     file, sampleName, sampleMass_mg, method, heatRate_Kmin (first segment),
%     segTable  struct with fields seg, descriptor (cellstr), rate_Kmin, gas
%     data      struct with column vectors seg, t_s, Tsample_C, Tref_C,
%               hf_mW, hfn_Wg
%
%   The header is parsed line by line; the numeric block is read in a single
%   textscan call, which is ~20x faster than readtable on these files.

fid = fopen(fpath, 'r', 'n', 'UTF-8');
if fid < 0, error('readJuliaDSC:open', 'Cannot open %s', fpath); end
c = onCleanup(@() fclose(fid));

S = struct('file', fpath, 'sampleName', '', 'sampleMass_mg', NaN, ...
           'method', '', 'heatRate_Kmin', NaN);
segs = {}; nHdr = 0; inSeg = false;
while true
    ln = fgetl(fid);
    if ~ischar(ln), error('readJuliaDSC:format', 'No data header in %s', fpath); end
    nHdr = nHdr + 1;
    ln = strtrim(strrep(ln, char(13), ''));
    if isempty(ln), continue, end
    if ln(1) ~= '#'
        if strncmpi(ln, 'Segment index', 13), break, end
        continue
    end
    ln = strtrim(ln(2:end));
    kv = regexp(ln, '^([^:\t]+):\s*(.*)$', 'tokens', 'once');
    if inSeg && ~isempty(regexp(ln, '^\d+\s', 'once'))
        segs(end+1, 1) = {ln}; %#ok<AGROW>
        continue
    end
    inSeg = false;
    if isempty(kv), continue, end
    key = lower(strtrim(kv{1})); val = strtrim(kv{2});
    if strncmp(key, 'sample name', 11),      S.sampleName = val;
    elseif strncmp(key, 'sample mass', 11),  S.sampleMass_mg = str2double(val);
    elseif strcmp(key, 'method'),            S.method = val;
    elseif strcmp(key, 'segments'),          inSeg = true;
    end
end

% numeric block: 7 numeric columns followed by a text column (protective gas)
C = textscan(fid, '%f %f %f %f %f %f %f %*s', 'Delimiter', '\t', ...
             'CollectOutput', true, 'ReturnOnError', true);
M = C{1};
M = M(all(isfinite(M(:,1:6)), 2), :);
S.data = struct('seg', M(:,1), 't_s', M(:,2), 'Tsample_C', M(:,3), ...
                'Tref_C', M(:,4), 'hf_mW', M(:,5), 'hfn_Wg', M(:,6));

% segment descriptors, e.g. "3  180degC/1(K/min)/300degC/Air"
n = numel(segs);
seg = nan(n,1); desc = cell(n,1); rate = nan(n,1); gas = cell(n,1);
for k = 1:n
    tok = regexp(segs{k}, '^(\d+)\s+(.*)$', 'tokens', 'once');
    seg(k) = str2double(tok{1}); desc{k} = strtrim(tok{2});
    r = regexp(desc{k}, '/([-\d.]+)\s*\(K/min\)', 'tokens', 'once');
    if ~isempty(r), rate(k) = str2double(r{1}); end
    p = regexp(desc{k}, '/', 'split');
    gas{k} = strtrim(p{end});
end
if n == 0                               % no descriptor block: infer from data
    seg = unique(S.data.seg); desc = repmat({''}, numel(seg), 1);
    rate = nan(numel(seg),1); gas = desc;
end
S.segTable = struct('seg', seg, 'descriptor', {desc}, 'rate_Kmin', rate, 'gas', {gas});
if ~isempty(rate), S.heatRate_Kmin = rate(1); end
end
