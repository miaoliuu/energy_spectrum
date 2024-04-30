clc; close all; clear;

%%

% 读取 速度场
data_uv = load("./vel_2d/vel_2d.mat");
u = data_uv.u; % unit[ m/s ]
v = data_uv.v; % unit[ m/s ]

lx = 0.18 * pi; % 速度场在 x 方向上的尺寸 unit[ m ]
ly = 0.18 * pi; % 速度场在 y 方向上的尺寸 unit[ m ]

%%

[k_st, knorm, wavenumbers, energy] = cal_spectrum_2d(u, v, lx, ly);

%%

% 读取 理论能谱数据 用于绘制图像
data_wnn = load("./vel_2d/wnn_2d.mat");
wnn = data_wnn.wnn;
wnn_spectrum = data_wnn.wnn_whichspec;
loglog(wnn, wnn_spectrum, "k-");
hold on;
% 绘制 计算得到的能谱
loglog(wavenumbers, energy, "bo--");
loglog([k_st, k_st], [1e-7, 1e-2], "r--")
hold off;
legend("theoretical", "calculated")
xlim([8, 10000])
ylim([1e-7, 1e-2])
xlabel("wave number (k)")
ylabel("E(k)")
grid on;
