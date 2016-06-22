%%
% STEP4: Related Analyse

%%
%--------Ԥ������---------------------------
%�������ݼ���ز�����ʼ������
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
    disp(['Topic_',num2str(i),':   ',char(WO_R(Index_P2_T_S(1,i)))]);  %���ÿ��topic��milestone service
end

WriteTopics( W_T_R , BETA , WO_R , 10 , 0.7 , 4 , 'results/1.1_TopicMilestoneService.txt' );    %���������ļ�
%%
% 1.2 Topic Temporal Strength
% D2: service-[topic-service]
disp('******************************************************');
disp('1.2 Cauculate Temporal Strength for each Service...');

disp('STEP1: Get Service Time...');
% ����ֻ����975�������ù�ϵ�ķ����ʵ����Ϣ
service_date = zeros(leftServiceNum,3);           %���зֱ�����ꡢ�¡���
sum_nodate = 0;
for i=1:leftServiceNum
    % ��ǰ�������ܵķ����б��б��ΪleftServiceSet(i)
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
minDate = [3000,1,1];       %��Сʱ��
maxDate = [1,1,1];          %���ʱ��
for i=1:leftServiceNum
    if (service_date(i,1)~=0)   %�ų��쳣����
        flag_max = 0;
        flag_min = 0;
        %�����ʱ��Ƚ�
        tem = service_date(i,:) - maxDate;
        if (tem(1)>0)       %���Ҫ��
            flag_max = 1;
        elseif (tem(1)==0)  %�����ͬ
            if (tem(2)>0)   %�·�Ҫ��
                flag_max = 1;
            elseif (tem(2)==0)  %�·���ͬ
                if (tem(3)>0)   %����Ҫ��
                    flag_max = 1;
                end
            end
        end
        %����Сʱ��Ƚ�
        tem = minDate - service_date(i,:);
        if (tem(1)>0)       %���Ҫ��
            flag_min = 1;
        elseif (tem(1)==0)  %�����ͬ
            if (tem(2)>0)   %�·�Ҫ��
                flag_min = 1;
            elseif (tem(2)==0)  %�·���ͬ
                if (tem(3)>0)   %����Ҫ��
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
disp(['ʱ��������Ϊ��λ����',num2str(dateLength),'��.']);

disp('STEP3: Draw Time Line for each Topic...');
TTS = zeros(T_R, dateLength+1);

%%
average_TTS = zeros(1,T_R);
for i=1:T_R         %��ÿ������ѭ��
    for j=1:K_R     %�Ե�ǰ�����µ�ǰK_R��service������ѭ��
        %topic i, service [serviceIndex]
        %topic-service distribution: P2_T_S
        serviceIndex = Index_P2_T_S(j,i);
        serviceTime = datenum(service_date(serviceIndex,1),service_date(serviceIndex,2),service_date(serviceIndex,3)) - datenum(minDate(1),minDate(2),minDate(3)) + 1;
        TTS(i, serviceTime) = TTS(i, serviceTime) + Sorted_P2_T_S(j,i); %1��7��������֮ǰ��д����
    end
%     for j=1:dateLength+1
%        TTS(i,j) = TTS(i,j)/sum(TTS(i,:)); 
%     end
    TTS(i,:) = TTS(i,:)/sum(TTS(i,:));              %1��7��������֮ǰ��д����
end

% ����ÿ�������ƽ��ʱ��
for i=1:T_R
   [~, temT] = find(TTS(i,:));
   sum_temT = sum(TTS(i,:));
   average_TTS(i) = sum(temT.*TTS(i,temT)/sum_temT); 
end

save 'results/1.2_TopicTemporalStrength.mat' TTS average_TTS;
%%
%ͼ��Ļ���
%plot(TTS(1,:))
% for i=1:35
%     DrawTTS(TTS,average_TTS,i,dateLength+1,360);
% end

%   0119֮ǰ
%   DrawTTS2(TTS,average_TTS,8,15,dateLength+1,360);
 
%   0119
 DrawTTS3(TTS,average_TTS,8,15,dateLength+1,360,Index_P2_T_S,Sorted_P2_T_S,service_date,minDate);
 
%%
% 1.3 Topic Keywords
disp('******************************************************');
disp('1.3 Topic Keywords Reconstribution:');
disp('load Descroption.mat for use...');

T2_W = zeros(T_R, length(WO_S));    %D2�е�topic-word����wordΪD1�е�word tokens

for i=1:T_R
    %ÿ��������ѭ��
    for j=1:K_R
        %ÿ�������µ�ǰK_R����ص�service��word��¼�����ۼ�
        %topic2 i, service Index_P2_T_S(j,i)
        %service_weight = Sorted_P2_T_S(j,i)
        serviceIndex = leftServiceSet(Index_P2_T_S(j,i)); 
        serviceWeight = Sorted_P2_T_S(j,i);
        T2_W(i,:) = serviceWeight*SW(serviceIndex,:);
    end
    %��һ��
    T2_W(i,:) = T2_W(i,:)/sum(T2_W(i,:));
end

%���ÿ��topicǰ30������
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

WriteTopicKeywords(T2_W,WO_S,T_R,topwordnum,'results/1.3_TopicKeywords.txt');   %���������ļ�

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

%���������ļ�
fid = fopen( 'results/2.1_TopicImportance.txt' , 'W' );
fprintf( fid , 'Top %d Topics are:\n',k);
for i=1:k
    fprintf( fid , 'Topic_%d:\t%10f\n',index(i),probtopic_R(index(i)));
end
fclose(fid);

