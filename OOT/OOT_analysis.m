function results = OOT_analysis(varargin)
%OOT_ANALYSIS  Oxidation onset temperature (OOT) of polyolefins by dynamic DSC.
%
%   OOT_analysis                   interactive: pick files, check labels, run
%   OOT_analysis('Files', {...})   analyse the given files (dialog still shown
%                                  for label checking unless 'Interactive',false)
%   R = OOT_analysis(..., Name, Value) also returns the results struct
%
%   Name-value options
%     'Files'          cellstr of CSV paths                 (default: ask)
%     'OutDir'         output folder          (default: <data folder>/OOT_results)
%     'Interactive'    show the run dialog         (default: desktop available)
%     'ExportFigures'  write PNG (300 dpi) + PDF (vector)          (default: true)
%     'PerRunFigures'  one tangent-construction figure per run     (default: true)
%
%   Reads JULIA / Anton Paar exports, isolates the oxidising ramp (last segment
%   run in air/O2) and applies the extrapolated-tangent onset construction on
%   the TEMPERATURE axis (ASTM E2009; dynamic variant of ISO 11357-6).
%
%   The sample label is taken from the file name: the text between the
%   dd/mm/yyyy date and "_OOT", e.g.
%     20261024_0451_Measurement_JULIA_DSC_23102026_HDPE_M1053_virgin_OOT.csv
%     -> HDPE_M1053_virgin
%   Runs sharing a label are replicates; the sample-level result is the
%   replicate mean with its half-range. Labels, replicate IDs and the polymer
%   (which sets the reference melting enthalpy) can be edited in the dialog.
%
%   Figures
%     Fig1  oxidising ramps, replicates overlaid (one line style per replicate)
%     Fig2  onset region, ordinate expanded
%     Fig3  tangent construction, one panel per run (+ one file per run)
%     Fig4  OOT by sample: replicates, mean and half-range (2+ runs)
%
%   Outputs in OutDir: OOT_summary.csv, OOT_by_sample.csv, OOT_QC.txt, fig/
%
%   TWO CONTROLS WORTH READING BEFORE THE NUMBERS
%
%   1. Specimen mass. OOT is sensitive to specimen mass because oxygen must
%      diffuse into the melt: a thicker specimen has a lower surface-to-volume
%      ratio, oxidises under greater diffusion control, and reports a higher
%      onset. Replicates whose masses differ materially are not replicates of
%      the onset, and the QC block flags the spread.
%
%   2. Mass-entry error. The tangent construction is invariant under a
%      constant scaling of the heat flow, so a wrong entered mass does not
%      bias the onset temperature, but it biases every amplitude. The script
%      therefore integrates the melting endotherm of the inert pre-heat and
%      reports the implied crystallinity: for the same material this must be
%      reproducible, and a value outside the range usual for the polymer
%      indicates the mass, not the specimen.
%
%   Requires readJuliaDSC.m, parseRunName.m, computeOOT.m, meltingCheck.m
%   (same folder). Base MATLAB only, no toolboxes. R2019b or later.

%% ---------------- options and configuration ----------------------------
ip = inputParser;
ip.addParameter('Files', {});
ip.addParameter('OutDir', '');
ip.addParameter('Interactive', usejava('desktop'));
ip.addParameter('ExportFigures', true);
ip.addParameter('PerRunFigures', true);
ip.parse(varargin{:});
opts = ip.Results;
if ischar(opts.Files), opts.Files = {opts.Files}; end

cfg = defaultConfig();

%% ---------------- file selection ---------------------------------------
files = opts.Files(:);
if isempty(files)
    if ~opts.Interactive
        error('OOT:noFiles', 'No files given. Pass ''Files'' or run interactively.');
    end
    files = pickFiles(cfg);
    if isempty(files), results = []; fprintf('No files selected.\n'); return, end
end
runs = loadRuns(files);
if isempty(runs), results = []; return, end
if isempty(opts.OutDir)
    opts.OutDir = fullfile(fileparts(runs(1).file), 'OOT_results');
end

%% ---------------- run dialog: confirm labels / replicates / polymer -----
if opts.Interactive
    [runs, opts] = runDialog(runs, opts, cfg);
    if isempty(runs), results = []; fprintf('Cancelled.\n'); return, end
end
runs = runs([runs.use]);
if isempty(runs), results = []; fprintf('No runs ticked.\n'); return, end

%% ---------------- analyse ----------------------------------------------
n = numel(runs);
for k = 1:n
    S = runs(k).S;
    seg = S.segTable;

    % oxidising segment: last one whose descriptor names air or oxygen
    isOx = ~cellfun(@isempty, regexpi(seg.descriptor, 'air|o2|oxygen', 'once'));
    if any(isOx)
        segOx = seg.seg(find(isOx, 1, 'last'));
    else
        segOx = max(S.data.seg);
        warning('OOT:noGasTag', 'No air/O2 tag in %s; using the last segment.', runs(k).name);
    end
    m  = S.data.seg == segOx;
    tt = S.data.t_s(m); TT = S.data.Tsample_C(m);
    R  = computeOOT(tt, TT, S.data.hfn_Wg(m), cfg.opt);

    % heating rate of the oxidising ramp itself (central 60 %)
    sp = R.TEnd_C - R.TStart_C;
    w  = TT >= R.TStart_C + 0.2*sp & TT <= R.TStart_C + 0.8*sp;
    R.rampOx_Kmin = NaN;
    if nnz(w) > 10, pOx = polyfit(tt(w)/60, TT(w), 1); R.rampOx_Kmin = pOx(1); end

    % melting endotherm of the inert pre-heat, as a control on the mass
    poly = cfg.polymers(strcmp({cfg.polymers.name}, runs(k).polymer));
    if isempty(poly), poly = cfg.polymers(end); end
    m1 = S.data.seg == min(seg.seg);
    Mk = meltingCheck(S.data.t_s(m1), S.data.Tsample_C(m1), S.data.hfn_Wg(m1), ...
                      struct('exoUp', cfg.opt.exoUp, 'dHref_Jg', poly.dHref_Jg));

    runs(k).R = R; runs(k).M = Mk; runs(k).XcRange = poly.XcRange;
end

%% ---------------- grouping ---------------------------------------------
labels = {runs.label}.';
[~, first, sIdx] = unique(labels);            % groups in order of first appearance
[first, ord] = sort(first); grp(ord) = 1:numel(ord);
uSample = labels(first); sIdx = reshape(grp(sIdx), [], 1);
nS = numel(uSample);
repOrd = zeros(n,1);
for g = 1:nS, ks = find(sIdx == g); repOrd(ks) = 1:numel(ks); end
for k = 1:n
    runs(k).color = cfg.palette(mod(sIdx(k)-1, size(cfg.palette,1)) + 1, :);
    runs(k).style = cfg.styles{mod(repOrd(k)-1, numel(cfg.styles)) + 1};
end

%% ---------------- tables -----------------------------------------------
R = [runs.R]; M = [runs.M];
OOT = [R.OOT_C].';
col = @(s, f) reshape([s.(f)], [], 1);
sumNames = {'sample','replicate','polymer','acquired','file','mass_mg', ...
            'rampInert_Kmin','rampOx_Kmin','T_start_C','T_final_C','OOT_C', ...
            'T_departure_C','status','T_inflect_C','maxSlope_WgK','T_peak_C', ...
            'peakHeight_Wg','dHox_Jg','dHox_open','baselineNoise_Wg', ...
            'Tm_peak_C','dHm_Jg','Xc_pct'};
fl = [R.flags];
sumCols = {labels, {runs.rep}.', {runs.polymer}.', {runs.acquired}.', {runs.name}.', ...
           col(runs,'mass'), arrayfun(@(r) r.S.heatRate_Kmin, runs(:)), ...
           col(R,'rampOx_Kmin'), col(R,'TStart_C'), col(R,'TEnd_C'), OOT, ...
           col(R,'TDeparture_C'), {R.status}.', col(R,'TInflect_C'), ...
           col(R,'slopeMax_WgK'), col(R,'TPeak_C'), col(R,'peakHeight_Wg'), ...
           col(R,'dHox_Jg'), double(col(fl,'enthalpyOpen')), col(R,'baseNoise_Wg'), ...
           col(M,'TPeak_C'), col(M,'dHm_Jg'), col(M,'Xc_pct')};

OOTmean = nan(nS,1); OOThalf = nan(nS,1); nRep = zeros(nS,1);
for g = 1:nS
    v = OOT(sIdx == g); v = v(~isnan(v));
    nRep(g) = numel(v);
    if ~isempty(v)
        OOTmean(g) = mean(v);
        OOThalf(g) = 0.5*(max(v) - min(v));         % half-range; n = 2 as a rule
    end
end

if exist(opts.OutDir, 'dir') ~= 7, mkdir(opts.OutDir); end
writeCsv(fullfile(opts.OutDir, 'OOT_summary.csv'), sumNames, sumCols);
writeCsv(fullfile(opts.OutDir, 'OOT_by_sample.csv'), ...
         {'sample','n','OOT_mean_C','OOT_halfRange_K'}, {uSample, nRep, OOTmean, OOThalf});

fprintf('\n=== OOT summary =============================================================\n');
fprintf('%-24s %-4s %7s %-11s %7s %8s %7s %6s %6s\n', 'Sample', 'Rep', 'OOT_C', 'status', ...
        'Tpk_C', 'dHox', 'm_mg', 'Xc_%', 'beta');
for k = 1:n
    fprintf('%-24s %-4s %7.1f %-11s %7.1f %8.0f %7.2f %6.1f %6.2f\n', runs(k).label, runs(k).rep, ...
            OOT(k), R(k).status, R(k).TPeak_C, R(k).dHox_Jg, runs(k).mass, M(k).Xc_pct, R(k).rampOx_Kmin);
end
fprintf('\nSample-level OOT (mean of replicates, half-range as spread)\n');
for g = 1:nS
    fprintf('  %-24s n = %d   OOT = %6.1f +/- %.1f degC\n', uSample{g}, nRep(g), OOTmean(g), OOThalf(g));
end

%% ---------------- QC block ---------------------------------------------
qcLines = qcReport(runs, sIdx, uSample, OOThalf, cfg);
fprintf('\nQC notes\n--------\n%s\n', strjoin(qcLines, newline));
fid = fopen(fullfile(opts.OutDir, 'OOT_QC.txt'), 'w', 'n', 'UTF-8');
if fid > 0, fprintf(fid, '%s\n', qcLines{:}); fclose(fid); end

%% ---------------- figures ----------------------------------------------
figDir = fullfile(opts.OutDir, 'fig');
rr = [R.rampOx_Kmin];
rateTxt = sprintf('%.3g K min^{-1} in air', median(rr(~isnan(rr))));
F = gobjects(0);

% Fig 1: oxidising ramps
F(end+1) = newFigure('OOT - oxidising ramps', [17 12]);
ax = axes(F(end)); hold(ax, 'on');
for k = 1:n
    r = runs(k).R;
    plot(ax, r.T, r.y, runs(k).style, 'Color', runs(k).color, 'DisplayName', runName(runs(k), true));
    if ~isnan(r.OOT_C)
        plot(ax, r.OOT_C, r.baseLevel_Wg, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', runs(k).color, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end
styleAxes(ax, cfg, 'Sample temperature (^{\circ}C)', 'Normalised heat flow (W g^{-1}), exo \uparrow', ...
          'Dynamic DSC in air: oxidising ramp', rateTxt);
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9);
xlim(ax, [min([R.TStart_C]) max([R.TEnd_C])]);
saveFig(F(end), figDir, 'Fig1_OOT_ramps', opts);

% Fig 2: onset region
ok = ~isnan(OOT);
if any(ok)
    H  = max([R(ok).peakHeight_Wg]);
    yb = [R(ok).baseLevel_Wg];
    xl = [min(OOT(ok)) - 20, max(OOT(ok)) + 12];
    yl = [min(yb) - 0.12*H, max(yb) + 0.55*H];
else
    xl = [min([R.TStart_C]) max([R.TEnd_C])]; yl = [-0.1 0.3];
end
F(end+1) = newFigure('OOT - onset region', [17 12]);
ax = axes(F(end)); hold(ax, 'on');
for k = 1:n
    r = runs(k).R;
    plot(ax, r.T, r.y, runs(k).style, 'Color', runs(k).color, 'DisplayName', runName(runs(k), true));
    if ~isnan(r.OOT_C)
        plot(ax, xl, polyval(r.baseCoef, xl), ':', 'Color', 0.55*[1 1 1], 'LineWidth', 0.8, ...
             'HandleVisibility', 'off');
        plot(ax, r.OOT_C, r.baseLevel_Wg, 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', runs(k).color, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    end
end
styleAxes(ax, cfg, 'Sample temperature (^{\circ}C)', 'Normalised heat flow (W g^{-1}), exo \uparrow', ...
          'Onset region, ordinate expanded', 'circles: extrapolated-tangent onset (OOT)');
legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9);
xlim(ax, xl); ylim(ax, yl);
saveFig(F(end), figDir, 'Fig2_onset_region', opts);

% Fig 3: tangent construction, one panel per run
nc = min(3, n); nr = ceil(n/nc);
F(end+1) = newFigure('OOT - tangent construction', [10*nc 8.5*nr + 1.5]);
tl = tiledlayout(F(end), nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:n
    ax = nexttile(tl);
    plotConstruction(ax, runs(k), cfg, false);
end
xlabel(tl, 'Sample temperature (^{\circ}C)', 'FontSize', cfg.LBL);
ylabel(tl, 'Normalised heat flow (W g^{-1}), exo \uparrow', 'FontSize', cfg.LBL);
title(tl, {'Extrapolated-tangent onset (ASTM E2009)', ['\\fontsize{10}\\color[rgb]{0.4 0.4 0.4}' ...
          'dashed grey: extrapolated baseline  |  black: tangent at max. slope  |  ' ...
          'dotted: first departure  |  shaded: baseline fit window']}, 'FontSize', cfg.TTL);
saveFig(F(end), figDir, 'Fig3_tangent_construction', opts);

% one stand-alone figure per run, named after the sample; when exported these
% are drawn off-screen and closed so that a large batch does not flood the desktop
if opts.PerRunFigures
    for k = 1:n
        f = newFigure(sprintf('OOT - %s %s', runs(k).label, runs(k).rep), [16 12.5]);
        if opts.ExportFigures, set(f, 'Visible', 'off'); end
        ax = axes(f, 'OuterPosition', [0 0.04 1 0.96]);       % room for the file-name footer
        plotConstruction(ax, runs(k), cfg, true);
        xlabel(ax, 'Sample temperature (^{\circ}C)', 'FontSize', cfg.LBL);
        ylabel(ax, 'Normalised heat flow (W g^{-1}), exo \uparrow', 'FontSize', cfg.LBL);
        annotation(f, 'textbox', [0 0 1 0.04], 'String', runs(k).name, 'Interpreter', 'none', ...
                   'EdgeColor', 'none', 'FontSize', 7, 'Color', 0.45*[1 1 1], ...
                   'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        saveFig(f, figDir, ['OOT_' safeName(runs(k).label) '_' safeName(runs(k).rep)], opts);
        if opts.ExportFigures, close(f); else, F(end+1) = f; end %#ok<AGROW>
    end
end

% Fig 4: OOT by sample (only informative with more than one run)
if n > 1
    F(end+1) = newFigure('OOT - by sample', [max(11, 3.2*nS + 7) 11]);
    ax = axes(F(end)); hold(ax, 'on');
    for g = 1:nS
        c  = cfg.palette(mod(g-1, size(cfg.palette,1)) + 1, :);
        ks = find(sIdx == g);
        if ~isnan(OOTmean(g))
            errorbar(ax, g, OOTmean(g), OOThalf(g), 'o', 'Color', c, 'MarkerFaceColor', c, ...
                     'MarkerSize', 9, 'LineWidth', 1.4, 'CapSize', 12);
            text(ax, g + 0.2, OOTmean(g), sprintf('%.1f \\pm %.1f', OOTmean(g), OOThalf(g)), ...
                 'FontSize', 9, 'Color', c);
        end
        jit = linspace(-0.13, 0.13, max(numel(ks), 2));
        for q = 1:numel(ks)
            k = ks(q);
            if isnan(OOT(k)), continue, end
            plot(ax, g + jit(q) - 0.15, OOT(k), 's', 'MarkerSize', 6, 'MarkerFaceColor', 'w', ...
                 'MarkerEdgeColor', c, 'LineWidth', 1.0);
            text(ax, g + jit(q) - 0.15, OOT(k), [runs(k).rep '  '], 'FontSize', 8, ...
                 'HorizontalAlignment', 'right');
        end
    end
    set(ax, 'XTick', 1:nS, 'XTickLabel', cellfun(@texEsc, uSample, 'UniformOutput', false), ...
            'XTickLabelRotation', 20 * (nS > 3));
    xlim(ax, [0.4 nS + 0.8]);
    if any(~isnan(OOT)), ylim(ax, [min(OOT) - 2, max(OOT) + 2]); end
    styleAxes(ax, cfg, 'Sample', 'OOT (^{\circ}C)', 'Oxidation onset temperature', ...
              'squares: replicates; filled: mean \pm half-range');
    saveFig(F(end), figDir, 'Fig4_OOT_by_sample', opts);
end

fprintf('\nResults written to %s\n', opts.OutDir);
if opts.ExportFigures, fprintf('Figures (PNG + PDF) in %s\n', figDir); end

if nargout
    results = struct('runs', runs, 'samples', struct('label', uSample, 'n', num2cell(nRep), ...
                     'OOT_mean_C', num2cell(OOTmean), 'OOT_halfRange_K', num2cell(OOThalf)), ...
                     'outDir', opts.OutDir, 'figures', F, 'qc', {qcLines});
end
end

%% ======================================================================
%%  configuration
%% ======================================================================
function cfg = defaultConfig()
% Fixed six-colour palette, one colour per sample; line style cycles per
% replicate so that a replicate pair reads as one material.
cfg.palette = [0.839 0.153 0.157;    % red
               0.122 0.467 0.706;    % blue
               0.173 0.627 0.173;    % green
               1.000 0.498 0.055;    % orange
               0.580 0.404 0.741;    % purple
               0.090 0.745 0.812];   % teal
cfg.styles  = {'-', '--', ':', '-.'};

% Reference enthalpy of the 100 % crystalline polymer and the crystallinity
% range treated as physically usual (QC only). Polymer is guessed from the
% first token of the label and can be changed in the dialog.
cfg.polymers = struct( ...
    'name',     {'PP',    'HDPE',  'MDPE',  'LDPE',  'LLDPE', 'PE',    'Other'}, ...
    'dHref_Jg', {207,     293,     293,     293,     293,     293,     NaN}, ...
    'XcRange',  {[30 58], [55 80], [45 62], [25 50], [25 55], [25 80], [0 100]});

cfg.opt = struct();
cfg.opt.exoUp      = true;     % oxidation exotherm positive in this export
cfg.opt.TExclude_K = 5.0;      % skip the first 5 K of the ramp (settling)
cfg.opt.minRise_Wg = 0.05;     % detection floor for accepting an exotherm

% Acceptance thresholds used only for the QC messages.
cfg.qc.massTol_pct = 15;       % replicate specimen masses should agree to this
cfg.qc.spreadTol_K = 2.0;      % replicate onset range worth remarking on

cfg.LBL = 13; cfg.TTL = 14;    % axis-label and title point sizes
end

%% ======================================================================
%%  file handling and run dialog
%% ======================================================================
function files = pickFiles(~)
start = getpref('OOT_analysis', 'lastDir', pwd);
if exist(start, 'dir') ~= 7, start = pwd; end
[f, p] = uigetfile({'*.csv', 'JULIA DSC export (*.csv)'; '*.*', 'All files'}, ...
                   'Select OOT measurement file(s)  -  Ctrl/Shift-click for several', ...
                   start, 'MultiSelect', 'on');
if isequal(f, 0), files = {}; return, end
setpref('OOT_analysis', 'lastDir', p);
files = fullfile(p, cellstr(f)).';
files = files(:);
end

function runs = loadRuns(files, runs)
% read each file once; label, replicate and polymer come from the file name
if nargin < 2
    runs = struct('file', {}, 'name', {}, 'S', {}, 'label', {}, 'suffix', {}, 'rep', {}, ...
                  'polymer', {}, 'mass', {}, 'acquired', {}, 'use', {});
end
for k = 1:numel(files)
    f = files{k};
    if any(strcmp({runs.file}, f)), continue, end
    try
        S = readJuliaDSC(f);
    catch err
        warning('OOT:read', 'Skipping %s: %s', f, err.message);
        continue
    end
    [~, nm, ext] = fileparts(f);
    I = parseRunName([nm ext], S.sampleName);
    runs(end+1) = struct('file', f, 'name', [nm ext], 'S', S, 'label', I.label, ...
                         'suffix', I.suffix, 'rep', '', 'polymer', I.polymer, ...
                         'mass', S.sampleMass_mg, 'acquired', I.acquired, 'use', true); %#ok<AGROW>
end
% chronological order, then R1, R2, ... within each label for unassigned runs
[~, o] = sort({runs.acquired}); runs = runs(o);
for k = 1:numel(runs)
    if isempty(runs(k).rep)
        same = strcmp({runs.label}, runs(k).label) & ~cellfun(@isempty, {runs.rep});
        runs(k).rep = sprintf('R%d', nnz(same) + 1);
    end
end
end

function [runs, opts] = runDialog(runs, opts, cfg)
% Table of runs: tick to include, edit label / replicate / polymer.
fig = figure('Name', 'OOT analysis - check runs', 'NumberTitle', 'off', 'MenuBar', 'none', ...
             'ToolBar', 'none', 'Color', 'w', 'Units', 'pixels', 'Position', [100 100 1040 480], ...
             'CloseRequestFcn', @(s, ~) finish(s, false));
movegui(fig, 'center');
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.02 0.91 0.96 0.06], ...
          'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'FontSize', 10, 'String', ...
          ['Sample labels are read from the file name (text between the date and "_OOT"). ' ...
           'Runs with the same label are replicates. Edit label, replicate or polymer if needed.']);
tbl = uitable(fig, 'Units', 'normalized', 'Position', [0.02 0.22 0.96 0.68], ...
              'ColumnName', {'Use', 'Sample label', 'Replicate', 'Polymer', 'Mass (mg)', 'Acquired', 'File'}, ...
              'ColumnFormat', {'logical', 'char', 'char', {cfg.polymers.name}, 'numeric', 'char', 'char'}, ...
              'ColumnEditable', [true true true true false false false], ...
              'ColumnWidth', {40, 210, 70, 75, 70, 120, 400}, 'RowName', [], 'FontSize', 10);
outTxt = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.16 0.12 0.64 0.06], ...
                   'String', opts.OutDir, 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.02 0.115 0.13 0.055], ...
          'String', 'Output folder:', 'BackgroundColor', 'w', 'HorizontalAlignment', 'right');
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.81 0.12 0.17 0.06], ...
          'String', 'Browse...', 'Callback', @browseOut);
cbFig = uicontrol(fig, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [0.02 0.03 0.22 0.06], ...
                  'String', 'Save figures (PNG + PDF)', 'Value', opts.ExportFigures, 'BackgroundColor', 'w');
cbRun = uicontrol(fig, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [0.25 0.03 0.22 0.06], ...
                  'String', 'One figure per run', 'Value', opts.PerRunFigures, 'BackgroundColor', 'w');
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.50 0.025 0.15 0.07], ...
          'String', 'Add files...', 'Callback', @addFiles);
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.67 0.025 0.14 0.07], ...
          'String', 'Cancel', 'Callback', @(~, ~) finish(fig, false));
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.83 0.025 0.15 0.07], ...
          'String', 'Run analysis', 'FontWeight', 'bold', 'Callback', @(~, ~) finish(fig, true));
setappdata(fig, 'runs', runs);
setappdata(fig, 'ok', false);
refresh();
uiwait(fig);

if ~ishandle(fig) || ~getappdata(fig, 'ok')
    runs = [];
else
    sync();
    runs = getappdata(fig, 'runs');
    opts.OutDir        = strtrim(get(outTxt, 'String'));
    opts.ExportFigures = logical(get(cbFig, 'Value'));
    opts.PerRunFigures = logical(get(cbRun, 'Value'));
end
if ishandle(fig), delete(fig); end

    function refresh()
        r = getappdata(fig, 'runs');
        set(tbl, 'Data', [{r.use}.', {r.label}.', {r.rep}.', {r.polymer}.', ...
                          {r.mass}.', {r.acquired}.', {r.name}.']);
    end
    function sync()
        r = getappdata(fig, 'runs'); d = get(tbl, 'Data');
        for i = 1:numel(r)
            r(i).use = logical(d{i,1});
            r(i).label = strtrim(d{i,2}); r(i).rep = strtrim(d{i,3}); r(i).polymer = d{i,4};
            if isempty(r(i).label), r(i).label = r(i).name; end
        end
        setappdata(fig, 'runs', r);
    end
    function addFiles(~, ~)
        sync();
        f = pickFiles(cfg);
        if isempty(f), return, end
        setappdata(fig, 'runs', loadRuns(f, getappdata(fig, 'runs')));
        refresh();
    end
    function browseOut(~, ~)
        d = uigetdir(get(outTxt, 'String'), 'Output folder');
        if ~isequal(d, 0), set(outTxt, 'String', d); end
    end
end

function finish(fig, ok)
setappdata(fig, 'ok', ok);
uiresume(fig);
if ~ok, set(fig, 'Visible', 'off'); end
end

%% ======================================================================
%%  QC
%% ======================================================================
function L = qcReport(runs, sIdx, uSample, OOThalf, cfg)
L = {};
for k = 1:numel(runs)
    r = runs(k).R; m = runs(k).M; fl = r.flags; msg = {};
    if fl.weakExotherm
        msg{end+1} = sprintf('no exotherm above %.2f W/g up to %.0f degC -> OOT censored', ...
                             cfg.opt.minRise_Wg, r.TEnd_C); %#ok<AGROW>
    end
    if fl.truncated,      msg{end+1} = 'inflection near the end of the ramp (exotherm truncated)'; end %#ok<AGROW>
    if fl.peakIncomplete, msg{end+1} = 'exotherm peak not reached before the end of the ramp'; end %#ok<AGROW>
    if fl.shortBaseline
        msg{end+1} = sprintf('baseline window only %.0f K (onset close to the start of the ramp)', ...
                             diff(r.baseWin_C)); %#ok<AGROW>
    end
    if fl.enthalpyOpen
        msg{end+1} = sprintf(['trace still %.2f W/g above baseline at %.0f degC: dHox = %.0f J/g ' ...
                              'is a LOWER BOUND, not a closed integral'], r.endRise_Wg, r.TEnd_C, r.dHox_Jg); %#ok<AGROW>
    end
    if ~isnan(r.TDeparture_C) && ~isnan(r.OOT_C) && (r.OOT_C - r.TDeparture_C) > 3
        msg{end+1} = sprintf('gradual pre-onset rise: departure at %.1f degC, tangent onset at %.1f degC', ...
                             r.TDeparture_C, r.OOT_C); %#ok<AGROW>
    end
    xr = runs(k).XcRange;
    if m.ok && ~isnan(m.Xc_pct) && (m.Xc_pct < xr(1) || m.Xc_pct > xr(2))
        msg{end+1} = sprintf(['melting check: Xc = %.1f %% is outside %.0f-%.0f %% usual for %s, which points ' ...
                              'to the ENTERED MASS (%.2f mg) rather than the specimen; onset unaffected, ' ...
                              'amplitudes are not'], m.Xc_pct, xr(1), xr(2), runs(k).polymer, runs(k).mass); %#ok<AGROW>
    end
    if m.ok && m.endothermOpen
        msg{end+1} = 'melting endotherm has not closed before the end of the inert ramp'; %#ok<AGROW>
    end
    if isempty(msg), msg = {'none'}; end
    L{end+1} = sprintf('  %-24s %-4s : %s', runs(k).label, runs(k).rep, msg{1}); %#ok<AGROW>
    for q = 2:numel(msg)
        L{end+1} = sprintf('  %-24s %-4s   %s', '', '', msg{q}); %#ok<AGROW>
    end
end

% replicate-level checks
fmt = @(v, f) strjoin(arrayfun(@(x) sprintf(f, x), v(:).', 'UniformOutput', false), '/');
for g = 1:numel(uSample)
    ks = find(sIdx == g);
    if numel(ks) < 2, continue, end
    reps = {runs(ks).rep};
    if numel(unique(reps)) < numel(reps)
        L{end+1} = sprintf('  ! %s: duplicate replicate IDs (%s).', uSample{g}, strjoin(reps, '/')); %#ok<AGROW>
    end
    mm = [runs(ks).mass];
    spreadPct = 100*(max(mm) - min(mm))/mean(mm);
    if spreadPct > cfg.qc.massTol_pct
        L{end+1} = sprintf(['  ! %s: replicate masses differ by %.0f %% (%s mg). OOT rises with specimen ' ...
                            'mass under oxygen-diffusion control, so this spread is confounded with the ' ...
                            'onset difference.'], uSample{g}, spreadPct, fmt(mm, '%.2f')); %#ok<AGROW>
    end
    oo = arrayfun(@(r) r.R.OOT_C, runs(ks));
    if ~isnan(OOThalf(g)) && 2*OOThalf(g) > cfg.qc.spreadTol_K
        L{end+1} = sprintf('  ! %s: replicate onsets differ by %.1f K (%s degC).', uSample{g}, ...
                           2*OOThalf(g), fmt(oo, '%.1f')); %#ok<AGROW>
    end
    xx = arrayfun(@(r) r.M.Xc_pct, runs(ks)); xx = xx(~isnan(xx));
    if numel(xx) > 1 && (max(xx) - min(xx)) > 5
        L{end+1} = sprintf(['  ! %s: normalised melting enthalpy differs between replicates (Xc %s %%). ' ...
                            'For one material this is a mass-entry discrepancy, not a crystallinity ' ...
                            'difference.'], uSample{g}, fmt(xx, '%.1f')); %#ok<AGROW>
    end
    dd = unique(cellfun(@(s) s(1:min(10, end)), {runs(ks).acquired}, 'UniformOutput', false));
    if numel(dd) > 1
        L{end+1} = sprintf(['  ! %s: replicates acquired on different dates (%s); confirm the ' ...
                            'calibration is unchanged.'], uSample{g}, strjoin(dd, ', ')); %#ok<AGROW>
    end
end
rr = arrayfun(@(r) r.R.rampOx_Kmin, runs); rr = rr(~isnan(rr));
if ~isempty(rr) && (max(rr) - min(rr)) > 0.05
    L{end+1} = '  ! oxidising heating rates are not identical across runs; OOT shifts with heating rate.';
end
end

%% ======================================================================
%%  plotting helpers
%% ======================================================================
function f = newFigure(name, sizeCm)
f = figure('Color', 'w', 'Name', name, 'NumberTitle', 'off', 'Units', 'centimeters', ...
           'Position', [2 2 sizeCm]);
set(f, 'DefaultAxesFontName', 'Helvetica', 'DefaultTextFontName', 'Helvetica', ...
       'DefaultAxesFontSize', 11, 'DefaultAxesLineWidth', 0.9, 'DefaultAxesBox', 'on', ...
       'DefaultAxesTickDir', 'out', 'DefaultLineLineWidth', 1.5);
end

function styleAxes(ax, cfg, xl, yl, ttl, sub)
xlabel(ax, xl, 'FontSize', cfg.LBL);
ylabel(ax, yl, 'FontSize', cfg.LBL);
if nargin > 5 && ~isempty(sub)
    ttl = {ttl, ['\fontsize{10}\color[rgb]{0.4 0.4 0.4}' sub]};
end
title(ax, ttl, 'FontSize', cfg.TTL, 'FontWeight', 'normal');
grid(ax, 'on'); set(ax, 'GridAlpha', 0.12);
end

function plotConstruction(ax, run, cfg, showLegend)
% trace, baseline, tangent, onset and departure for one run, fully labelled
hold(ax, 'on');
r = run.R; c = run.color; grey = 0.45*[1 1 1];
if ~isnan(r.OOT_C)
    xr   = [r.OOT_C - 25, max(r.OOT_C + 10, min(r.TEnd_C, r.TPeak_C + 15))];
    span = max(0.25, 2.5*(r.yInflect_Wg - r.baseLevel_Wg));
    yl   = [r.baseLevel_Wg - 0.15*span, r.baseLevel_Wg + span];
else
    xr = [r.TStart_C r.TEnd_C];
    bl = polyval(r.baseCoef, r.T);
    yl = [min(bl) - 0.05, max(min(bl) + 0.25, max(r.y) + 0.05)];
end
h = gobjects(0);
h(end+1) = patch(ax, [r.baseWin_C fliplr(r.baseWin_C)], yl([1 1 2 2]), 0.93*[1 1 1], ...
                 'EdgeColor', 'none', 'DisplayName', 'Baseline fit window');
h(end+1) = plot(ax, r.T, r.y, '-', 'Color', c, 'LineWidth', 1.7, 'DisplayName', 'Heat flow (air)');
h(end+1) = plot(ax, xr, polyval(r.baseCoef, xr), '--', 'Color', grey, 'LineWidth', 1.1, ...
                'DisplayName', 'Extrapolated baseline');
if ~isnan(r.OOT_C)
    tx = [r.OOT_C - 6, r.TInflect_C + (yl(2) - r.yInflect_Wg)/r.slopeMax_WgK];
    h(end+1) = plot(ax, tx, r.yInflect_Wg + r.slopeMax_WgK*(tx - r.TInflect_C), '-', ...
                    'Color', [0.15 0.15 0.15], 'LineWidth', 1.1, 'DisplayName', 'Tangent at max. slope');
    if ~isnan(r.TDeparture_C)
        h(end+1) = plot(ax, r.TDeparture_C*[1 1], yl, ':', 'Color', grey, 'LineWidth', 1.1, ...
                        'DisplayName', 'First departure');
    end
    plot(ax, r.OOT_C*[1 1], yl, '--', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    h(end+1) = plot(ax, r.OOT_C, r.baseLevel_Wg, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w', ...
                    'MarkerEdgeColor', c, 'LineWidth', 1.6, 'DisplayName', 'OOT (onset)');
    info = {sprintf('\\bfOOT = %.1f ^{\\circ}C\\rm', r.OOT_C), ...
            sprintf('T_{dep} = %.1f ^{\\circ}C', r.TDeparture_C), ...
            sprintf('T_{peak} = %.1f ^{\\circ}C', r.TPeak_C)};
    if ~strcmp(r.status, 'ok'), info{end+1} = sprintf('status: %s', r.status); end
else
    info = {sprintf('\\bfOOT > %.0f ^{\\circ}C\\rm', r.TEnd_C), 'no exotherm (censored)'};
end
text(ax, 0.03, 0.96, info, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 0.8*[1 1 1], 'Margin', 4);
xlim(ax, xr); ylim(ax, yl);
sub = sprintf('%s  |  %.2f mg  |  %.3g K min^{-1} in air', run.acquired, run.mass, r.rampOx_Kmin);
title(ax, {['\bf' texEsc(run.label) '\rm  \cdot  ' texEsc(run.rep)], ...
           ['\fontsize{9}\color[rgb]{0.4 0.4 0.4}' sub]}, 'FontSize', 12, 'FontWeight', 'normal');
grid(ax, 'on'); set(ax, 'GridAlpha', 0.12, 'Layer', 'top');
if showLegend
    legend(ax, h, 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 3, ...
           'FontSize', 9, 'Box', 'off');
end
end

function s = runName(run, withOOT)
s = sprintf('%s  \\cdot  %s', texEsc(run.label), texEsc(run.rep));
if withOOT
    if isnan(run.R.OOT_C)
        s = sprintf('%s   (OOT > %.0f ^{\\circ}C)', s, run.R.TEnd_C);
    else
        s = sprintf('%s   (OOT %.1f ^{\\circ}C)', s, run.R.OOT_C);
    end
end
end

function s = texEsc(s)
% escape TeX specials so that e.g. HDPE_M1053_virgin is shown literally
s = regexprep(s, '([_^{}\\])', '\\$1');
end

function s = safeName(s)
s = regexprep(s, '[^\w\-]', '_');
end

function saveFig(f, figDir, name, opts)
drawnow;
if ~opts.ExportFigures, return, end
if exist(figDir, 'dir') ~= 7, mkdir(figDir); end
base = fullfile(figDir, name);
if exist('exportgraphics', 'file') || exist('exportgraphics', 'builtin')
    exportgraphics(f, [base '.png'], 'Resolution', 300);
    exportgraphics(f, [base '.pdf'], 'ContentType', 'vector');
else
    set(f, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'auto');
    pp = get(f, 'PaperPosition'); set(f, 'PaperSize', pp(3:4));
    print(f, [base '.png'], '-dpng', '-r300');
    print(f, [base '.pdf'], '-dpdf', '-painters');
end
end

%% ======================================================================
%%  CSV writer (no table dependency, UTF-8, quoted text)
%% ======================================================================
function writeCsv(fname, names, cols)
n = numel(cols{1});
fid = fopen(fname, 'w', 'n', 'UTF-8');
if fid < 0, warning('OOT:write', 'Cannot write %s', fname); return, end
fprintf(fid, '%s\n', strjoin(names, ','));
for i = 1:n
    cells = cell(1, numel(cols));
    for j = 1:numel(cols)
        v = cols{j};
        if iscell(v)
            cells{j} = ['"' strrep(v{i}, '"', '""') '"'];
        elseif isnan(v(i))
            cells{j} = '';
        else
            cells{j} = sprintf('%.6g', v(i));
        end
    end
    fprintf(fid, '%s\n', strjoin(cells, ','));
end
fclose(fid);
end
