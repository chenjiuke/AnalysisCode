function info = parseRunName(fname, headerName)
%PARSERUNNAME  Extract the sample label and acquisition time from a file name.
%
%   info = parseRunName('20261024_0451_Measurement_JULIA_DSC_23102026_HDPE_M1053_virgin_OOT.csv')
%     info.label    = 'HDPE_M1053_virgin'   text between the dd/mm/yyyy date and "_OOT"
%     info.suffix   = ''                    text after "_OOT" (e.g. 'replicate')
%     info.acquired = '2026-10-24 04:51'    from the leading yyyymmdd_HHMM stamp
%     info.dateKey  = '20261024'
%     info.polymer  = 'HDPE'                inferred from the first token of the label
%
%   If the file name does not follow the pattern, the "#Sample name" header
%   (second argument, optional) is tried, then the bare file name is used.

if nargin < 2, headerName = ''; end
[~, base] = fileparts(fname);
info = struct('label', '', 'suffix', '', 'acquired', '', 'dateKey', '', 'polymer', 'Other');

% leading acquisition stamp yyyymmdd_HHMM
tok = regexp(base, '^(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})', 'tokens', 'once');
if ~isempty(tok)
    info.dateKey  = [tok{1:3}];
    info.acquired = sprintf('%s-%s-%s %s:%s', tok{:});
end

% <anything>_<8-digit date>_<LABEL>_OOT[_<suffix>]
pat = '^(?:.*_)?\d{8}_(.+?)_OOT(?:_(.+))?$';
tok = regexp(base, pat, 'tokens', 'once', 'ignorecase');
if isempty(tok) && ~isempty(headerName)
    tok = regexp(headerName, pat, 'tokens', 'once', 'ignorecase');
end
if isempty(tok)
    % last resort: drop the stamp and instrument prefix, keep the rest
    lab = regexprep(base, '^\d{8}_\d{4}_(Measurement_)?(JULIA_DSC_)?', '');
    tok = {regexprep(lab, '_OOT$', '', 'ignorecase'), ''};
end
info.label = tok{1};
if numel(tok) > 1 && ~isempty(tok{2}), info.suffix = tok{2}; end

% polymer family from the first token of the label
first = upper(regexp(info.label, '^[A-Za-z]+', 'match', 'once'));
% (exact match only, so e.g. 'PET' is not taken for PE; editable in the dialog)
if any(strcmp(first, {'PP','PE','HDPE','MDPE','LDPE','LLDPE'}))
    info.polymer = first;
elseif strcmp(first, 'UHMWPE')
    info.polymer = 'HDPE';
end
end
