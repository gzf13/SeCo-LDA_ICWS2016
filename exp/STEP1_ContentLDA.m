%%
%experiment1
%edit by gzf 1209
%STEP1:Content-based LDA for Service


%%
%�������ݼ�������ʼ��

load('../Description.mat');

%%
%��ȡÿһ�����ʵ���Ϣ�����浽�ĵ���
fid = fopen('step1_wordcount_service.txt','wt');
for i=1:1:sNum
    %Service_i
    temIndex = find(SW(i,:));
    temL = length(temIndex);
    for j=1:1:temL
        %������ļ� 
        fprintf(fid,'%d %d %d\n',i,temIndex(j),SW(i,temIndex(j)));
    end
end
fclose(fid);

%%
% Content-based LDA
% ����׼��
[WS1, DS1] = importworddoccounts('step1_wordcount_service.txt');    %����importworddoccounts�����������
WO_S = WD';                                                   %ת��1*W�ĵ����ֵ����

% Set the number of topics
T=35; 

% Set the hyperparameters
BETA=0.01;
ALPHA=50/T;

% The number of iterations
N = 400; 

% The random seed
SEED = 3;

% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 2;

% This function might need a few minutes to finish
tic
%WP DP Z
%DIC: WO_S
[ W_T_S, S_T, Z_S ] = GibbsSamplerLDA( WS1 , DS1 , T , N , ALPHA , BETA , SEED , OUTPUT );
toc

%%
%Ԥ������һ������D1�������
%DIC: WO_S
%[WP,DP,Z]
%[ W_T_S, S_T, Z_S ] = GibbsSamplerLDA( WS1 , DS1 , T , N , ALPHA , BETA , SEED , OUTPUT );
%D1��service-topic-word

%ÿ������ĸ���
W_S = size( W_T_S , 1 );
T_S = size( W_T_S , 2 );
sumW_T_S = sum( W_T_S , 1 ) + BETA*W_S;
probtopic = sumW_T_S / sum( sumW_T_S ); 

%ÿ�������µ��ʵĸ���
%ǰ100��
K_S = 100;                              %��ʾǰ100������
Sorted_P1_T_W = zeros( K_S , T_S );     %D1��topic-word���ʷֲ��ľ���ֵ��С���Ӵ�С��
Index_P1_T_W = zeros( K_S , T_S );      %D1��topic-word���ʷֲ���word��ǩ
%���е���
for t=1:T_S
   [ temp1 , temp2 ] = sort( -W_T_S( : , t ) );
   Sorted_P1_T_W( : , t )  = ( full( -temp1( 1:K_S )) + BETA ) ./ ( repmat( sumW_T_S( t ) , K_S , 1 ));
   Index_P1_T_W( : , t )   = temp2( 1:K_S );
end

P1_T_W = zeros( length(WO_S) , T_S );             %D1��ÿ��������word�ķֲ������D1_topic_wordֵ�Ĵ�С��δ����,���п�
for t=1:T_S
    temp1 = W_T_S(:,t);
    P1_T_W( : , t )  = ( temp1 + BETA ) / sumW_T_S( t );
end

P1_S_T = zeros(sNum, T_S);                %D1��ÿ���ĵ�(service)������ֲ������D1_service_topicֵ�ô�С�����п�
for i=1:sNum
    P1_S_T(i,:) = (S_T(i,:) + ALPHA)/(sum(S_T(i,:)) + ALPHA*T_S);
    i
end

%%
% �������������ļ�
save STEP1DATA Sorted_P1_T_W Index_P1_T_W P1_T_W P1_S_T WO_S sNum probtopic



