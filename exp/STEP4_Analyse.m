%%
% STEP4: Related Analyse

%%
%--------预备工作---------------------------
%载入数据及相关参数初始化设置
clear;clc;
load 'STEP1DATA.mat';
load 'STEP2DATA.mat';
load 'STEP3DATA.mat';
load '../Record.mat';
load '../Description.mat'
T_S = 35;
T_R = 35;
BETA=0.01;
%%
% -----------------------------------------------------------------
% PART1: Single Topic
% -----------------------------------------------------------------

%%
% 1.1 Topic Milestone Service
% D2: service-[topic-service]
disp('******************************************************');
disp('1.1 Topic Milestone Service:');
for i=1:T_R
    disp(['Topic_',num2str(i),':   ',char(WO_R(Index_P2_T_S(1,i)))]);  %输出每个topic的milestone service
end

WriteTopics( W_T_R , BETA , WO_R , 10 , 0.7 , 4 , 'results/1.1_TopicMilestoneService.txt' );    %保存结果到文件
%%
% 1.2 Topic Temporal Strength
% D2: service-[topic-service]
disp('******************************************************');
disp('1.2 Cauculate Temporal Strength for each Service...');

disp('STEP1: Get Service Time...');
% 这里只计算975个有引用关系的服务的实践信息
service_date = zeros(leftServiceNum,3);           %三列分别代表年、月、日
sum_nodate = 0;
for i=1:leftServiceNum
    % 当前服务在总的服务列表中编号为leftServiceSet(i)
    tem = service(leftServiceSet(i),4);
    if (isempty(tem{1})==1)
        sum_nodate = sum_nodate+1;
    else
        temS = char(service(leftServiceSet(i),4));
        dotIndex = strfind(temS,'.');
        temL = length(temS);
        service_date(i,1) = str2num(temS((dotIndex(2)+1):temL));
        service_date(i,2) = str2num(temS(1:(dotIndex(1)-1)));
        service_date(i,3) = str2num(temS((dotIndex(1)+1):(dotIndex(2)-1)));
    end
end

%%
disp('STEP2: Definite axis...');
minDate = [3000,1,1];       %最小时间
maxDate = [1,1,1];          %最大时间
for i=1:leftServiceNum
    if (service_date(i,1)~=0)   %排除异常数据
        flag_max = 0;
        flag_min = 0;
        %与最大时间比较
        tem = service_date(i,:) - maxDate;
        if (tem(1)>0)       %年份要大
            flag_max = 1;
        elseif (tem(1)==0)  %年份相同
            if (tem(2)>0)   %月份要大
                flag_max = 1;
            elseif (tem(2)==0)  %月份相同
                if (tem(3)>0)   %天数要大
                    flag_max = 1;
                end
            end
        end
        %与最小时间比较
        tem = minDate - service_date(i,:);
        if (tem(1)>0)       %年份要大
            flag_min = 1;
        elseif (tem(1)==0)  %年份相同
            if (tem(2)>0)   %月份要大
                flag_min = 1;
            elseif (tem(2)==0)  %月份相同
                if (tem(3)>0)   %天数要大
                    flag_min = 1;
                end
            end
        end
        if (flag_max==1)    
            maxDate = service_date(i,:); 
        end
        if (flag_min==1)    
            minDate = service_date(i,:);
        end
    end
end
disp(['MaxDate:',num2str(maxDate)]);
disp(['MinDate:',num2str(minDate)]);
dateLength = datenum(maxDate(1),maxDate(2),maxDate(3)) - datenum(minDate(1),minDate(2),minDate(3));
disp(['时间轴以天为单位，共',num2str(dateLength),'天.']);

disp('STEP3: Draw Time Line for each Topic...');
TTS = zeros(T_R, dateLength+1);

%%
average_TTS = zeros(1,T_R);
for i=1:T_R         %对每个主题循环
    for j=1:K_R     %对当前主题下的前K_R个service依次做循环
        %topic i, service [serviceIndex]
        %topic-service distribution: P2_T_S
        serviceIndex = Index_P2_T_S(j,i);
        serviceTime = datenum(service_date(serviceIndex,1),service_date(serviceIndex,2),service_date(serviceIndex,3)) - datenum(minDate(1),minDate(2),minDate(3)) + 1;
        TTS(i, serviceTime) = TTS(i, serviceTime) + Sorted_P2_T_S(j,i); %1月7号修正，之前编写错误。
    end
%     for j=1:dateLength+1
%        TTS(i,j) = TTS(i,j)/sum(TTS(i,:)); 
%     end
    TTS(i,:) = TTS(i,:)/sum(TTS(i,:));              %1月7号修正，之前编写错误。
end

% 计算每个主题的平均时间
for i=1:T_R
   [~, temT] = find(TTS(i,:));
   sum_temT = sum(TTS(i,:));
   average_TTS(i) = sum(temT.*TTS(i,temT)/sum_temT); 
end

save 'results/1.2_TopicTemporalStrength.mat' TTS average_TTS;
%%
%图像的绘制
%plot(TTS(1,:))
% for i=1:35
%     DrawTTS(TTS,average_TTS,i,dateLength+1,360);
% end

%   0119之前
%   DrawTTS2(TTS,average_TTS,8,15,dateLength+1,360);
 
%   0119
 DrawTTS3(TTS,average_TTS,8,15,dateLength+1,360,Index_P2_T_S,Sorted_P2_T_S,service_date,minDate);
 
%%
% 1.3 Topic Keywords
disp('******************************************************');
disp('1.3 Topic Keywords Reconstribution:');
disp('load Descroption.mat for use...');

T2_W = zeros(T_R, length(WO_S));    %D2中的topic-word矩阵，word为D1中的word tokens

for i=1:T_R
    %每个主题做循环
    for j=1:K_R
        %每个主题下的前K_R个相关的service的word记录依次累加
        %topic2 i, service Index_P2_T_S(j,i)
        %service_weight = Sorted_P2_T_S(j,i)
        serviceIndex = leftServiceSet(Index_P2_T_S(j,i)); 
        serviceWeight = Sorted_P2_T_S(j,i);
        T2_W(i,:) = serviceWeight*SW(serviceIndex,:);
    end
    %归一化
    T2_W(i,:) = T2_W(i,:)/sum(T2_W(i,:));
end

%输出每个topic前30个单词
topwordnum = 50;
% for i=1:T_R
%     %topic i in D2
%     [value, indexOrder] = sort(-T2_W(i,:));
%     disp('-----------------------------------------------------------');
%     disp(['Topic_',num2str(i),' Keywords:',num2str(sum(T2_W(i,indexOrder(1:topwordnum))))]);
%     for j=1:topwordnum
%        disp([char(WO_S(indexOrder(j))),'     ',num2str(T2_W(i,indexOrder(j)))]); 
%     end
% end

WriteTopicKeywords(T2_W,WO_S,T_R,topwordnum,'results/1.3_TopicKeywords.txt');   %输出结果到文件

%%
%-----------------------------------
% PART2:Service System

%%
% 2.1 Topic Importance
disp('******************************************************');
disp('2.1 Topic Importance:');

[value,index] = sort(-probtopic_R);
k=35;
disp(['Top ',num2str(k),' Topics are:']);
for i=1:k
   disp(['Topic_',num2str(index(i)),':     ',num2str(probtopic_R(index(i)))]); 
end

%输出结果到文件
fid = fopen( 'results/2.1_TopicImportance.txt' , 'W' );
fprintf( fid , 'Top %d Topics are:\n',k);
for i=1:k
    fprintf( fid , 'Topic_%d:\t%10f\n',index(i),probtopic_R(index(i)));
end
fclose(fid);

