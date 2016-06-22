%%
% edit on 12.21
% STEP5 make recommendation for service-service usage

%%
%%
%--------预备工作---------------------------
%载入数据及相关参数初始化设置
load 'STEP1DATA.mat';
load 'STEP2DATA.mat';
load 'STEP3DATA.mat';
load '../Record.mat';
load '../Description.mat'
BETA = 0.01;
ALPHA = 50/35;

%%
% Service-Service Citation结果
% 求出每两个service之间的使用相关度，即P(i,j)=Pr(si->sj|si);
ss_use = zeros(leftServiceNum,leftServiceNum);      %对应的数据下标存在leftServiceSet中
for i = 1:leftServiceNum
    %求出每个service的service使用相关度
    % 已知i，使用j
    p_si_t = P2_S_T(i,:);%1*50
    p_t_sj = P2_T_S;   % sNum_R*50 
    ss_use(i,:) = p_t_sj * p_si_t';  %相乘结果为leftServiceNum*1; i - other services
    i
end
%ss_use归一化
for i=1:leftServiceNum
    ss_use(i,:) = ss_use(i,:)/sum(ss_use(i,:));
end

%%
% Mashup-Service Citation 结果
load 'STEP5_0_DATA.mat';
ss_use2 = zeros(msSNum_new,msSNum_new);             %对应的数据下标存在ssum中

% TIPS: P(M=j|S=i)=1/sum(msRelation(i,:));

%%
% 求出每两个service之间的使用相关度，即P(i,j)=Pr(si->sj|si);
for i = 1:msSNum_new
    %求出每个service的service使用相关度
    % 已知i，使用j
    for j=1:msSNum_new
        % P(sj|si)=sigma[k]P(sj|topic=k)*P(topic=k|si)
        sum_p_sj_si = 0;
        for k=1:T_MSR
           % topic k
           % P(sj|topic=k)
           p_sj_tk = P3_T_S(j,k);
           
           % P(toipc=k|si) = sigma[mashup=n]P(topic=k|mashup=n)*P(mashup=n|si)
           p_tk_mn = P3_M_T(:,k);   %mNum*1
           p_mn_si = msRelation(i,:)/sum(msRelation(i,:));  %1*mNum
           p_tk_si = sum(p_tk_mn.*p_mn_si');
           
           sum_p_sj_si = sum_p_sj_si + p_sj_tk*p_tk_si;
        end
        ss_use2(i,j) = sum_p_sj_si;
    end
    i
end

%ss_use归一化
for i=1:msSNum_new
    ss_use2(i,:) = ss_use2(i,:)/sum(ss_use2(i,:));
end

%%
% *********************************************************************************************************************
% 验证一：对引用关系的还原程度
% *********************************************************************************************************************
%%
% part1: service-service还原service-service
threshold = 0.00001:0.00001:0.005;
deltaSum = zeros(1,length(threshold));
ssRelationOri = ssRelationLeft./repmat(sum(ssRelationLeft,2),1,leftServiceNum); %归一化原始的引用矩阵
for i = 1:length(threshold)
    flagMartrix = ss_use>=threshold(i);
    ssRelationRB = ss_use.*flagMartrix;
    % 至少保证一个服务有一个推荐结果，否则将误差输出为-1
    if (sum(sum(ssRelationRB,2)==0)>0)
        deltaSum(i) = -1;
    else
        % 归一化
        ssRelationRB = ssRelationRB./repmat(sum(ssRelationRB,2),1,leftServiceNum);
        % 计算误差
        deltaSum(i) = sum(sum((ssRelationRB - ssRelationOri).*(ssRelationRB - ssRelationOri)));
        deltaSum(i) = deltaSum(i)/leftServiceNum;
    end
    threshold(i)
end
[value,index] = min(deltaSum);
disp('********service-service还原***************')
disp(['最小误差：',num2str(value)]);
disp(['相应阈值：',num2str(threshold(index))]);
disp('******************************************');
figure();plot(threshold,deltaSum);
bestT = threshold(index);
%%
% part2: mashup-service还原service-service
threshold2 = 0.00011:0.00001:0.001;
deltaSum2 = zeros(1,length(threshold2));
for i = 1:length(threshold2)
    flagMartrix = ss_use2>=threshold2(i);
    ssRelationRB = ss_use2.*flagMartrix;
    % 至少保证一个服务有一个推荐结果，否则将误差输出为-1
    if (sum(sum(ssRelationRB,2)==0)>0)
        deltaSum2(i) = -1;
    else
        % 归一化
        ssRelationRB = ssRelationRB./repmat(sum(ssRelationRB,2),1,leftServiceNum);
        % 计算误差
        deltaSum2(i) = sum(sum((ssRelationRB - ssRelationOri).*(ssRelationRB - ssRelationOri)));
        deltaSum2(i) = deltaSum2(i)/leftServiceNum;
    end
    threshold2(i)
end
[value,index] = min(deltaSum2);
disp('*********mashup-service还原*************')
disp(['最小误差：',num2str(value)]);
disp(['相应阈值：',num2str(threshold2(index))]);
disp('***************************************');
figure();plot(threshold2,deltaSum2);
bestT2 = threshold2(index);


%%
% *********************************************************************************************************************
% 验证二：单个Service推荐结果MAP@N结果测试
% *********************************************************************************************************************
%%
% part1: service-service还原service-service
% part2: mashup-service还原service-service
% part3: Apriori还原service-service
load('STEP6_Apriori_0217.mat');
% part4: service content还原service-service
load('STEP6_CMSD_0217.mat');

% 根据相应阈值得到最佳的矩阵
flagMartrix = ss_use>=bestT;
ssRelationRB = ss_use.*flagMartrix;

flagMartrix2 = ss_use2>=bestT2;
ssRelationRB2 = ss_use2.*flagMartrix2;

% 对原始记录ssRelationOri进行排序
Sorted_ssRelationOri = zeros(leftServiceNum,leftServiceNum);
for i=1:leftServiceNum
    [~, Sorted_ssRelationOri(i,:)] = sort(-ssRelationOri(i,:));
end
% 对service-service还原结果ssRelationRB进行排序
Sorted_ssRelationRB = zeros(leftServiceNum,leftServiceNum);
for i=1:leftServiceNum
    [~, Sorted_ssRelationRB(i,:)] = sort(-ssRelationRB(i,:));
end
% 对mashup-service还原结果ssRelationRB2进行排序
Sorted_ssRelationRB2 = zeros(leftServiceNum,leftServiceNum);
for i=1:leftServiceNum
    [~, Sorted_ssRelationRB2(i,:)] = sort(-ssRelationRB2(i,:));
end
% 对Apriori还原结果ssRelationRB3进行排序
Sorted_ssRelationRB3 = zeros(leftServiceNum,leftServiceNum);
for i=1:leftServiceNum
    [~, Sorted_ssRelationRB3(i,:)] = sort(-p_l2(i,:));
end
% 对service content还原结果ssRelationRB3进行排序
Sorted_ssRelationRB4 = zeros(leftServiceNum,leftServiceNum);
for i=1:leftServiceNum
    [~, Sorted_ssRelationRB4(i,:)] = sort(-P4(i,:));
end


% 计算MAP@N
N = 1:1:10;
result_map = zeros(1,length(N));
result_map2 = zeros(1,length(N));
result_map3 = zeros(1,length(N));
result_map4 = zeros(1,length(N));
for i=1:length(N)
    % MAP@ni
    ni = N(i);
    for j=1:leftServiceNum
       result_map(i) = result_map(i) + calMAP_N(Sorted_ssRelationOri(j,:),Sorted_ssRelationRB(j,:),ni);
       result_map2(i) = result_map2(i) + calMAP_N(Sorted_ssRelationOri(j,:),Sorted_ssRelationRB2(j,:),ni); 
       result_map3(i) = result_map3(i) + calMAP_N(Sorted_ssRelationOri(j,:),Sorted_ssRelationRB3(j,:),ni); 
       result_map4(i) = result_map4(i) + calMAP_N(Sorted_ssRelationOri(j,:),Sorted_ssRelationRB4(j,:),ni); 
    end
    result_map(i) = result_map(i)/leftServiceNum;
    result_map2(i) = result_map2(i)/leftServiceNum;
    result_map3(i) = result_map3(i)/leftServiceNum;
    result_map4(i) = result_map4(i)/leftServiceNum;
    N(i)
end
plot(result_map,'--rs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','r');
hold on
plot(result_map3,'--gs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','g');
plot(result_map4,'--ms','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','m');
plot(result_map2,'--bs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','b');
legend('SeCo-LDA','AA','CMSD','MUR-LDA');
xlabel('N');
ylabel('MAP');
grid on;
%%
% *********************************************************************************************************************
% 验证三：平均最小JS距离对比、困惑度对比
% *********************************************************************************************************************
load 'STEP5_0_scitation.mat';
minD_JS_SCitation = minD_JS;
per_SCitation = perplexity;
load 'STEP5_0_mcitation.mat';
minD_JS_MCitation = minD_JS;
per_MCitation = perplexity;
load 'STEP5_0_scontext.mat';
minD_JS_SContext = minD_JS;
per_SContext = perplexity;
T_N = 5:5:70;
plot(T_N,minD_JS_SCitation,'--bs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','b');
hold on
plot(T_N,minD_JS_MCitation,'--rs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','r');
plot(T_N,minD_JS_SContext,'--gs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','g');
title('AVERAGE MIN D_J_S');
legend('service citation','masup citation','service context');
%%
figure();
plot(T_N,per_SCitation,'--bs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','b');
hold on
plot(T_N,per_MCitation,'--rs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','r');
plot(T_N,per_SContext,'--gs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','g');
title('perplexity');
legend('service citation','masup citation','service context');
