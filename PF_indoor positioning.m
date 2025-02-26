clear all;
M = 200;
K = 300;
% M = 100;
% K = 100;

%% Set up of geometry location
nn = 3; AA = -52;
% rssi = xlsread('Data.xlsx','final_med');
rssi = xlsread('Data.xlsx','final_coordinate');
d = 10.^((abs(rssi(:,3:5)-AA)/(10*nn)))*1;
d = d(:,[3 1 2]);


a = 0;   b = 2*pi;
N = 3;
Speed_of_Light = 3e8;
sample_rate = 61.44e8;
Tc = 1/sample_rate;
A = (-1*diag(ones(1,N-1))+diag(ones(1,N-2),1));
A = [A [zeros(N-2,1);1]];

L_length = 0.925;
W_length = 0.775;

% 
gNB_Position =[-3*L_length 0*L_length 3*L_length  ;
               0           6*W_length 0*W_length  ];
theta = angle(gNB_Position(1,:)+gNB_Position(2,:)*1j);

% Exponential form
gNB_Position_exp = gNB_Position(1,:)+1j*gNB_Position(2,:);


%% First step: Initialization sample particles from initial distribution
% We establish the initial point of 2D-plane
% Initial x_predict
x_est = zeros(4,1);

mse_LS = [];
mse_PF = [];
 
%% Particle implement
x_Est_record = [];
x_Real_record  = [rssi(:,1)*L_length   + -3*L_length ...
                rssi(:,2)*W_length   + 0.00000000001   ]';
x_Real_record_exp  = x_Real_record(1,:)+1j*x_Real_record(2,:);
for ID = 1:25
    UE_Position  = [rssi(ID,1)*L_length   + -3*L_length;...
                    rssi(ID,2)*W_length   + 0.00000000001   ];
    UE_Position_exp  = UE_Position(1)+1j*UE_Position(2);
    x_P = zeros(4,1) + 1*randn(4,M);
    q_ = ones(1,M)/M;
    for k =1:K
    %     UE_Position = tracking_line(:,mod(k,80)+1);
    %     UE_Position_exp  = UE_Position(1)+1j*UE_Position(2);
    
        % ToA & TDoA & DoA
        DoA = d(ID,:)';
        ToA = DoA/Speed_of_Light;
        
        %% Secd step: Uppdate & Prediction
        % Update
        x_P_update = x_P + 0.1*randn(size(x_P));
        Dist = [];
        for m=1:M
            Dist = [Dist sqrt(sum(((x_P_update(1:2,m)-gNB_Position)).^2,1))'];
        end
    
        %% Third step: Update evaluate the pdf with measurement
        for m=1:M
    %         q(m) = 1/sum(abs(DoA-Dist(:,m)).^2);
            q(m) = q_(m)*exp(-sum(abs(DoA-Dist(:,m)).^2)/(2*20));
        end
        
        %% Normaliza the relative likelihood
        % Sort order
        q_ = q/sum(q);
    
        %% Fourth step: Resample
        for m=1:M
            Index_of_resampled = find(rand <= cumsum(q_),1);
            x_P(:,m) = x_P_update(:,Index_of_resampled) + 0.0005*randn(4,1);
        end
    
%         x_est = mean(x_P,2);
        x_est = sum(q_.*x_P,2);
        x_est_exp = x_est(1)+1j*x_est(2);
    
        mse_PF = [mse_PF mean((x_est(1:2)-UE_Position).^2)];
    
        
%         %% Plot geometry location
%         figure(1)
%         plot(x_P(1,:)+x_P(2,:)*1j,'o');
%         hold on;
%         grid on;
%         plot(x_est(1,:)+x_est(2,:)*1j,'ro','LineWidth',3);
%         plot(gNB_Position_exp,'k+','LineWidth',3);
%         plot(UE_Position(1)+UE_Position(2)*1j,'bo','LineWidth',4);
%         plot(x_Real_record(1,:)+x_Real_record(2,:)*1j,'b*')
%     
%         %% ToA circle plot
%         for iii = 1:N
%             angles = linspace(0, 2*pi, 500);
%             radius = mean(ToA(iii,:))*Speed_of_Light;
%             CenterX = gNB_Position(1,iii);
%             CenterY = gNB_Position(2,iii);
%             x = radius * cos(angles) + CenterX;
%             y = radius * sin(angles) + CenterY;
%             plot(CenterX, CenterY, 'LineWidth', 3, 'MarkerSize', 14);
%             plot(x, y, 'LineWidth', 2);
%         end
%     
%         legend('Particle','PF Est','gNB','UE');
%         title('Geometry Location');
%         xlim([-10 10])
%         ylim([-5 10])
%         hold off;
        
    end
    x_Est_record = [x_Est_record x_est(1:2)];
end

hold on;
MSE = sqrt(sum((x_Est_record-x_Real_record).^2,1))
mean(MSE)
plot(MSE)
title('Error distance v.s. ID')
legend(['Reference power A is ',num2str(AA)]);
xlabel('ID of Sampled Point');
ylabel('Error distance in meter');