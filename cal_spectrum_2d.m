function [knyquist, knorm, kappa_kr, tke_kr] = cal_spectrum_2d(u, v, lx, ly)
    % 输入参数：u, v - 速度场在 x, y 方向上的分量
    % 输入参数：lx, ly - 速度场在 x, y 方向上的尺寸
    % 返回值：knyquist - 有效波数截止区间
    % 返回值：knorm - 正则化波数
    % 返回值：kappa_kr - 球形积分后的波数
    % 返回值：tke_kr - 球形积分后的能量
    if not( all(size(u) == size(v)) )
        error("u and v needs same size!")
    end
    [nx, ny] = size(u);
    nt = nx * ny;
    % step 1 对速度场进行傅里叶变换
    uf = fftshift(fft2(u)) / nt;
    vf = fftshift(fft2(v)) / nt;
    % step 1.1 傅里叶变换后对应的波数空间 (wave number, kappa)
    kappa_x = floor(-nx/2) : 1 : floor(nx/2)-1;
    kappa_y = floor(-ny/2) : 1 : floor(ny/2)-1;
    % step 2 速度场能量谱 (kinetic energy spectrum of velocity-sepctral)
    tkef = 0.5 .* real(uf .* conj(uf) + vf .* conj(vf));
    % step 3 通过在波数空间的球面积分进行归一化 (surface integral at kappa space)
    tke_kr = zeros([1, max([nx, ny])]);
    for i = 1:nx
        for j = 1:ny
            kx = kappa_x(i);
            ky = kappa_y(j);
            kr = round(sqrt(kx .* kx + ky .* ky));
            tke_kr(kr+1) = tke_kr(kr+1) + tkef(i, j);
        end
    end
    k0x = 2.0 .* pi ./ lx;
    k0y = 2.0 .* pi ./ ly;
    knorm = (k0x + k0y) ./ 2.0;
    knyquist = knorm .* min([nx, ny]) ./ 2;
    kappa_kr = knorm .* (0:1:length(tke_kr)-1);
    tke_kr = tke_kr ./ knorm;
end
