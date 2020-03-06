% 基于Matlab的二指夹具在手操作仿真
clear;  clc;

% --- 符号规定: 
% g: 重力加速度, t: 时间序列, h: 位置序列, hd: 目标位置, 
% hv: 位置误差, dhv: 速度误差, s: 跟踪误差,

% --- 控制参数调整
lamb = 1;
k = 10;
alpha_a = 4;
alpha_b = 4;
alpha_c = 4;
ae(1) = 0;
be(1) = 0;
ce(1) = 0;

% --- 模型常量估计值
g = 9.8;
miu = 0.5;
sigma = 0.5;
m = 1;

% --- 系统初值及参考模型
tt = 0.01;
t = 0: tt: 10;
N = length(t);
h0 = 0;
hd = 1+ones(N,1);

% --- 循环变量初始化
h = zeros(N,1);
dh = zeros(N,1);
ddh = zeros(N,1);

% For begining
h(1) = h0;
dh(1) = 0;
ddh(1) = 0;
dhv(1) = 0;

% Loop begin
% 1:N-1 is better
for i = 2:N
    h(i) =  h(i-1) + dh(i-1)*tt + 1/2*ddh(i-1)*tt*tt;
    
    dh(i) = h(i) - h(i-1);
    dhd(i) = hd(i) - hd(i-1);
    ddhd(i) = dhd(i) - dhd(i-1);
    hr(i) = ddhd(i) - lamb*dhv;
    hv = h(i) - hd(i);
    dhv = dh(i) - dhd(i);
    ddhv = ddh(i) - ddhd(i);
    s(i) = dhv + lamb*hv;
    
    ae(i) = ae(i-1) - alpha_a*s(i)*hr(i-1);
    be(i) = be(i-1) - alpha_b*s(i)*dh(i-1);
    ce(i) = ce(i-1) - alpha_c*s(i)*g;
    
    ufn(i) = ae(i)*hr(i) - k*s(i) +be(i)*dh(i) + ce(i)*g;
    fn(i) = ufn(i);
    ddh(i) = miu/m*fn(i) - sigma/m*dh(i) - g;
end

plot(t, hd, ':');
hold on;
plot(t, h);

legend('Ref', 'Real');
