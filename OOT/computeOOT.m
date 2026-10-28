function R = computeOOT(t, T, y, opt)
%COMPUTEOOT  Extrapolated-tangent oxidation onset on the temperature axis.
%
%   R = computeOOT(t_s, T_C, y_Wg, opt) analyses one oxidising ramp
%   (ASTM E2009, dynamic variant of ISO 11357-6):
%
%     1. the ramp is binned onto a uniform temperature grid (opt.dT_K), which
%        averages out sampling noise and makes all later windows K-based;
%     2. a straight baseline is fitted to the pre-onset region, iterating the
%        window so that it ends opt.baseGap_K before the first departure;
%     3. the tangent is drawn at the steepest point of the leading edge;
%     4. OOT = intersection of the baseline and that tangent.
%
%   opt fields (all optional)
%     exoUp       true    oxidation exotherm is positive in the export
%     TExclude_K  5       skip this much of the start of the ramp (settling)
%     minRise_Wg  0.05    exotherm must exceed baseline by this to be accepted
%     depThr_Wg   0.01    floor of the "first departure" threshold
%     nSigma      5       departure threshold in units of baseline noise
%     dT_K        0.1     temperature grid step
%     smooth_K    1.0     moving-mean width applied before differentiating
%     deriv_K     2.0     moving-mean width applied to dy/dT
%     baseLen_K   30      maximum baseline window length
%     baseGap_K   3       gap between baseline window and departure point
%
%   R contains the onset quantities, the binned trace (R.T, R.y, R.t) and the
%   construction (baseCoef, baseWin_C, TInflect_C, yInflect_Wg, slopeMax_WgK)
%   needed to redraw it, plus R.flags and a one-word R.status.

def = struct('exoUp', true, 'TExclude_K', 5, 'minRise_Wg', 0.05, 'depThr_Wg', 0.01, ...
             'nSigma', 5, 'dT_K', 0.1, 'smooth_K', 1.0, 'deriv_K', 2.0, ...
             'baseLen_K', 30, 'baseGap_K', 3);
if nargin < 4, opt = struct(); end
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opt, fn{k}), opt.(fn{k}) = def.(fn{k}); end
end

t = t(:); T = T(:); y = y(:);
if ~opt.exoUp, y = -y; end

% ---------- bin onto a uniform temperature grid ------------------------------
T0 = min(T) + opt.TExclude_K;  T1 = max(T);
keep = T >= T0;  t = t(keep); T = T(keep); y = y(keep);
dT = opt.dT_K;
ib  = floor((T - T0)/dT) + 1;
cnt = accumarray(ib, 1);
ok  = cnt > 0;
Tb  = accumarray(ib, T) ./ max(cnt,1);
yb  = accumarray(ib, y) ./ max(cnt,1);
tb  = accumarray(ib, t) ./ max(cnt,1);
Tg  = (T0 + dT/2 : dT : T1).';
yg  = interp1(Tb(ok), yb(ok), Tg, 'linear', 'extrap');
tg  = interp1(Tb(ok), tb(ok), Tg, 'linear', 'extrap');
n   = numel(Tg);

ys = movmean(yg, max(1, round(opt.smooth_K/dT)));
ds = movmean(gradient(ys, dT), max(1, round(opt.deriv_K/dT)));

% ---------- iterate baseline window / departure point ------------------------
win = [Tg(1), min(Tg(1) + 20, Tg(1) + (Tg(end) - Tg(1))/3)];
for it = 1:5
    w = Tg >= win(1) & Tg <= win(2);
    winUsed = win;
    p = polyfit(Tg(w), yg(w), 1);
    r = ys - polyval(p, Tg);
    noise = std(yg(w) - polyval(p, Tg(w)));
    [pk, ip] = max(r);
    if pk < opt.minRise_Wg, break, end
    [~, im] = max(ds(1:ip));
    thr = max(opt.nSigma*noise, opt.depThr_Wg);
    j = find(r(1:im) <= thr, 1, 'last');
    if isempty(j), j = 1; end
    newWin = [max(Tg(1), Tg(j) - opt.baseGap_K - opt.baseLen_K), Tg(j) - opt.baseGap_K];
    if diff(newWin) < 5                         % onset close to ramp start
        newWin = [Tg(1), max(Tg(1) + 2, Tg(j) - 0.5)];
    end
    if all(abs(newWin - win) < dT), break, end
    win = newWin;
end

R = struct('T', Tg, 'y', yg, 't', tg, 'TStart_C', Tg(1), 'TEnd_C', Tg(end), ...
           'baseCoef', p, 'baseWin_C', winUsed, 'baseNoise_Wg', noise, ...
           'OOT_C', NaN, 'TDeparture_C', NaN, 'TInflect_C', NaN, 'yInflect_Wg', NaN, ...
           'slopeMax_WgK', NaN, 'baseLevel_Wg', NaN, 'TPeak_C', NaN, ...
           'peakHeight_Wg', pk, 'dHox_Jg', NaN, 'endRise_Wg', r(end), 'status', 'ok');
R.flags = struct('weakExotherm', pk < opt.minRise_Wg, 'truncated', false, ...
                 'peakIncomplete', false, 'shortBaseline', diff(winUsed) < 10, ...
                 'enthalpyOpen', false);
if R.flags.weakExotherm
    R.status = 'censored';
    return
end

% ---------- tangent construction ---------------------------------------------
s  = ds(im);
yi = ys(im);
Ti = Tg(im);
if s <= p(1)                                    % tangent never meets baseline
    R.status = 'censored'; R.flags.weakExotherm = true;
    return
end
R.OOT_C        = (yi - s*Ti - p(2)) / (p(1) - s);
R.TDeparture_C = Tg(j);
R.TInflect_C   = Ti;
R.yInflect_Wg  = yi;
R.slopeMax_WgK = s;
R.baseLevel_Wg = polyval(p, R.OOT_C);
R.TPeak_C      = Tg(ip);

% ---------- oxidation enthalpy: departure -> return to baseline (or ramp end)
k2 = ip - 1 + find(r(ip:end) <= thr, 1, 'first');
if isempty(k2), k2 = n; end
R.dHox_Jg = trapz(tg(j:k2), r(j:k2));           % W/g * s = J/g

R.flags.truncated      = Ti > Tg(end) - 2;
R.flags.peakIncomplete = ip > n - round(1/dT);
R.flags.enthalpyOpen   = k2 == n && r(end) > max(thr, 0.05*pk);
if R.flags.truncated
    R.status = 'truncated';
elseif R.flags.peakIncomplete
    R.status = 'provisional';
end
end
