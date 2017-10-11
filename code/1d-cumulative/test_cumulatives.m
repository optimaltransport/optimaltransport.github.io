%%
% Test for 1D OT using cumulative distributions.

addpath('../toolbox/');

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
SetTickOff = @()set(gca, 'XTick',[], 'YTick',[]);
SetTickOn = @()set(gca, 'XTick',[0 1/2 1], 'YTick',[0 1/2 1]);

rep = 'results/1d-cumulative/';
[~,~] = mkdir(rep);

N = 256;
t = linspace(0,1,N)';
gauss = @(m,s)exp(-(t-m).^2/(2*s).^2);
normalize = @(x)x/sum(x(:));

s = .6;
a = gauss(.3*s,.05*s) + .5*gauss(.6*s,.15*s);
b = .5*gauss(1-s*.2,.04*s) + .8*gauss(1-s*.5,.05*s) + .5*gauss(1-s*.8,.04*s);
vmin = .025;
a = normalize(vmin + a);
b = normalize(vmin + b);

plot(t,[a b]);

clf;
area(t, a, 'FaceColor', 'r');
axis tight; SetAR(1/2); SetTickOff();
saveas(gca, [rep 'input-mu.eps'], 'epsc');

clf;
area(t, b, 'FaceColor', 'b');
axis tight; SetAR(1/2); SetTickOff();
saveas(gca, [rep 'input-nu.eps'], 'epsc');

% cumulative
ca = cumsum(a);
cb = cumsum(b);
% inverse cumulatives
ica = interp1(ca, t, t, 'spline');
icb = interp1(cb, t, t, 'spline');
% composition of function
Tab = interp1(t, icb, ca, 'spline'); % icb o ca
Tba = interp1(t, ica, cb, 'spline'); % ica o cb
% should be close to Id
Iaa = interp1(t, Tba, Tab, 'spline'); % Tba o Tab
Ibb = interp1(t, Tab, Tba, 'spline'); % Tab o Tba

lw = 2;
clf; hold on;
plot(t, ca, 'r', 'LineWidth', lw);
plot(t, cb, 'b', 'LineWidth', lw);
axis tight; box on; SetTickOn();  SetAR(1);
legend('C_{\mu}', 'C_{\nu}');
set(gca, 'FontSize', 20);
saveas(gca, [rep 'cumul.eps'], 'epsc');

clf; hold on;
plot(t, ica, 'r', 'LineWidth', lw);
plot(t, icb, 'b', 'LineWidth', lw);
axis tight; box on; SetTickOn();  SetAR(1);
legend('C_{\mu}^{-1}', 'C_{\nu}^{-1}');
set(gca, 'FontSize', 20);
saveas(gca, [rep 'icumul.eps'], 'epsc');

clf; hold on;
plot(t, Tab, 'color', [1/2 1/2 0], 'LineWidth', lw);
plot(t, Tba, 'color', [0 1/2 1/2], 'LineWidth', lw);
legend('T', 'T^{-1}');
axis tight; box on; SetTickOn();  SetAR(1);
set(gca, 'FontSize', 20);
saveas(gca, [rep 'transports.eps'], 'epsc');

% Barycenter of inverse cumulant`
q = 9;
m = {}; icm = {}; cm = {};
for i=1:q
    r = (i-1)/(q-1);
    icm{i} = (1-r)*ica + r*icb;
    cm{i} = interp1(icm{i}, t, t, 'spline');
    m{i} = diff([0;cm{i}]);
end

% display inverpolation of inverse cumulant
clf; hold on;
for i=1:q
    r = (i-1)/(q-1);
    col = [1-r 0 r];
    plot(t, icm{i}, 'color', col, 'LineWidth', lw);
end
axis tight; SetAR(1); SetTickOn(); box on;
set(gca, 'FontSize', 20);
saveas(gca, [rep 'interp-cumul.eps'], 'epsc');

% display inverpolation of measures
clf; hold on;
for i=1:q
    r = (i-1)/(q-1);
    col = [1-r 0 r];
    plot(t, m{i}, 'color', col, 'LineWidth', lw);
end
axis tight; SetAR(1/2); SetTickOff(); box on;
set(gca, 'FontSize', 20);
saveas(gca, [rep 'interp-bary.eps'], 'epsc');


