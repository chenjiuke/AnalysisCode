function M = meltingCheck(t, T, y, opt)
%MELTINGCHECK  Integrate the melting endotherm of the inert pre-heat.
%
%   M = meltingCheck(t_s, T_C, y_Wg, opt) returns dHm (J/g), the implied
%   crystallinity Xc = dHm/dHref, and the melting peak temperature. Because
%   the OOT construction is invariant to a mis-entered specimen mass while
%   every amplitude is not, a reproducible Xc is the check on the mass entry.
%
%   opt fields (all optional)
%     exoUp     true          sign convention of the export
%     dHref_Jg  NaN           enthalpy of the 100 % crystalline polymer
%     intWin_C  []            [Tlow Thigh] integration limits; [] = automatic
%                             (55 K below to 18 K above the endotherm peak,
%                             clipped to the ramp)
%     skip_K    10            ignore the start-up transient of the ramp
%
%   The baseline is the straight line joining the trace (2 K means) at the
%   two integration limits. M.endothermOpen is set when the upper limit is
%   clipped by the end of the ramp and the trace there is still changing.

if nargin < 4, opt = struct(); end
if ~isfield(opt, 'exoUp'),    opt.exoUp = true;   end
if ~isfield(opt, 'dHref_Jg'), opt.dHref_Jg = NaN; end
if ~isfield(opt, 'intWin_C'), opt.intWin_C = [];  end
if ~isfield(opt, 'skip_K'),   opt.skip_K = 10;    end

M = struct('ok', false, 'dHm_Jg', NaN, 'Xc_pct', NaN, 'TPeak_C', NaN, ...
           'intWin_C', [NaN NaN], 'endothermOpen', false);
t = t(:); T = T(:); y = y(:);
if numel(T) < 50, return, end
if opt.exoUp, y = -y; end                   % endotherm positive from here on

use = T >= min(T) + opt.skip_K;
t = t(use); T = T(use); y = y(use);
ys = movmean(y, max(1, round(numel(y)/500)));
Tlo = min(T); Thi = max(T);

[~, ip] = max(ys);
M.TPeak_C = T(ip);
if isempty(opt.intWin_C)
    win = [max(Tlo + 2, M.TPeak_C - 55), min(Thi - 1, M.TPeak_C + 18)];
    clipped = M.TPeak_C + 18 > Thi - 1;
else
    win = [max(Tlo, opt.intWin_C(1)), min(Thi, opt.intWin_C(2))];
    clipped = opt.intWin_C(2) > Thi;
end
if diff(win) < 10, return, end
M.intWin_C = win;

a = abs(T - win(1)) <= 1;  b = abs(T - win(2)) <= 1;
y1 = mean(ys(a)); y2 = mean(ys(b));
w  = T >= win(1) & T <= win(2);
bl = y1 + (y2 - y1) * (T(w) - win(1)) / diff(win);
M.dHm_Jg = trapz(t(w), ys(w) - bl);        % W/g * s = J/g
if isfinite(opt.dHref_Jg) && opt.dHref_Jg > 0
    M.Xc_pct = 100 * M.dHm_Jg / opt.dHref_Jg;
end

% still changing at the upper limit? compare end slope with pre-melt slope
if clipped
    e = T >= win(2) - 3 & T <= win(2);
    s = T >= win(1) & T <= win(1) + 10;
    pe = polyfit(T(e), ys(e), 1); ps = polyfit(T(s), ys(s), 1);
    M.endothermOpen = abs(pe(1)) > 3*abs(ps(1)) + 0.005;
end
M.ok = true;
end
