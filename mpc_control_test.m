%**************************************************************************
 % @file       		mpc_control_test.m
 % @company         SZU
 % @author     		黎国彬
 % @Software 		Matlab2016a
 % @Target         `利用MPC跟踪简单轨迹测试
 % @date            2022-3-1
 % All rights reserved
%**************************************************************************

clear;clc;
close all;

%% 选择是否生成gif
gif_generate_flag = 0;                                                     %若要将结果导出为gif，则将该flag置1
pic_num = 1;

%% 仿真参数设定
dt = 0.01;                                                                 % 离散时间间隔
num_step = 1500;                                                           % 仿真步数
t = (0:num_step-1)*dt;                                                     % 生成仿真时间序列

%% 运动学参数设定
L = 1;                                                                     %轴距
%设置参考输入(参考速度vel_ref和参考前轮转角delta_ref)
%u_ref = [3;0]*ones(1,num_step);                                           %相应参考轨迹为直线
u_ref = [1;0.5]*ones(1,num_step);                                          %相应参考轨迹为圆                                        

%设定参考状态，即参考轨迹，state_ref的三行分别存储x_ref，y_ref和theta_ref 
state_ref = zeros(3,num_step);                                              
for k = 2:num_step
    last_theta = state_ref(3,k-1);
    state_ref(:,k) = state_ref(:,k-1)+[cos(last_theta),0;sin(last_theta),0;0,1]*u_ref(:,k-1)*dt;
end

%设定初始状态和初始输入
state_real = zeros(3,num_step);                                            %初始化实际状态向量
state_real(:,1) = [1;-1;0];                                                %设定初始状态
state_error = zeros(3,num_step);                                           %初始化状态误差向量
u_real = zeros(2,num_step);                                                %初始化实际输入向量
u_real(:,1) = [0;0];                                                       %设定初始输入
u_error = zeros(2,num_step);                                               %初始化输入误差向量

%% MPC控制器参数设定
N = 10;                                                                    %预测区间
Q = diag([1,1,0.5]);                                                       %状态量权重
R = diag([0.1,0.1]);                                                       %控制量权重
%计算Q_bar和R_bar（权重矩阵）
Q_bar = kron(eye(N+1),Q);
R_bar = kron(eye(N),R);

%% 二次型优化器设置（使用 interior-point-convex 算法，不显示迭代过程）
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off');

%% 滚动优化
tic;                                                                       %开始计时
figure(1);
hold on;
xlabel('x(m)');
ylabel('y(m)');
%为了不处理最后N个时刻的特殊情况，循环只进行到 num_step-N
for k = 1:num_step-N
   %刻画每个时刻的参考位置和实时位置
   plot(state_real(1,k),state_real(2,k),'b.');
   plot(state_ref(1,k),state_ref(2,k),'r.');
   drawnow;                                                                %实时显示轨迹
   %录制gif动图
   if gif_generate_flag
       F=getframe(gcf);
        I=frame2im(F);
        [I,map]=rgb2ind(I,256);
        if pic_num == 1
            imwrite(I,map,'MPC_test.gif','gif', 'Loopcount',inf,'DelayTime',0.00);
        else
            imwrite(I,map,'MPC_test.gif','gif','WriteMode','append','DelayTime',0.00);
        end
        pic_num = pic_num + 1;
   end
   %计算当前时刻的状态误差
   state_error(:,k) = state_real(:,k) - state_ref(:,k);
   %计算A和B
   theta_ref = state_ref(3,k);
   A = eye(3)+dt*[0,0,-u_ref(1,k)*sin(theta_ref);0,0,u_ref(1,k)*cos(theta_ref);0,0,0];
   B = dt*[cos(theta_ref),0;sin(theta_ref),0;tan(u_ref(2,k))/L,u_ref(1,k)*(sec(u_ref(2,k)))^2/L];
   n = size(A,1);                                                          %A是nxn矩阵，求出n
   p = size(B,2);                                                          %B为nxp矩阵，求出p
   %计算M和C
   M = [eye(n);zeros(N*n,n)];                                              %初始化M矩阵
   C = zeros((N+1)*n,N*p);                                                 %初始化C矩阵
   tmp = eye(n);
   for i=1:N
       rows=i*n+(1:n);                                                     %定义当前行数，从i*n行开始，共n行
       C(rows,:) = [tmp*B,C(rows-n,1:end-p)];                              %填装C矩阵，后半部分代码是将前n行的内容复制并删接到当前n行的后部
       tmp = A*tmp;
       M(rows,:) = tmp;                                                    %填装M矩阵
   end
   %计算G,E,H
   G = M'*Q_bar*M;
   E = M'*Q_bar*C;
   H = C'*Q_bar*C+R_bar;
   %调用 quadprog 求解器进行最优化
   f = ((state_error(:,k))'*E)';
   [U_k,fval,exitflag,output,lambda] = quadprog(H,f,[],[],[],[],[],[],[],options);
   u_error(:,k)= U_k(1:2);                                                 %只取U_k的第一部分（前两项）得到输入误差
   u_real(:,k) = u_ref(:,k) + u_error(:,k);                                %实际输入 = 参考输入 + 误差输入
   state_error(:,k+1) = A*state_error(:,k)+B*u_error(:,k);                 %由当前时刻的状态误差和输入误差可求出下一时刻的预测状态误差
   state_real(:,k+1) = state_ref(:,k+1) + state_error(:,k+1);              %下一时刻状态 = 下一时刻参考状态 + 预测误差状态
end
toc;                                                                       %结束计时
