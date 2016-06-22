%%
%experiment3
%edit by gzf 1219
%STEP3:Relation-based LDA

%%
% Content-based LDA
% ����׼��

load('../Record.mat');
load('STEP2DATA.mat');
[WS, DS] = importworddoccounts('step2_wordcount_ssRelationLeft.txt');    %����importworddoccounts�����������

%ʹ����ϴ��������

%WOΪservice�б�
%ѡȡ�������ù�ϵ��service��Ϊ�ʵ�
%��Щservice���ܵ�service�µı�ż�¼��leftServiceSet�У�����ΪleftServiceNum

%WO = service(:,1);
WO = service(leftServiceSet,1);
WO_R = WO';

% Set the number of topics
T=35; 

% Set the hyperparameters
BETA=0.01;
ALPHA=50/T;

% The number of iterations
N = 500; 

% The random seed
SEED = 3;

% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 2;

% This function might need a few minutes to finish
tic
%[WP,DP,Z]
[ W_T_R,S_T_R,Z_R ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
toc

%%
%Ԥ��������������D2�������
%DIC: WO_R
%[WP,DP,Z]
%[ W_T_R,S_T_R,Z_R ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
%D2:service-topic-service

%ÿ������ĸ���
W_R = size( W_T_R , 1 );                    %���ʵ�������D2��Ϊservice������
T_R = size( W_T_R , 2 );                    %�����������50
sumW_T_R = sum( W_T_R , 1 ) + BETA*W_R;
probtopic_R = sumW_T_R / sum( sumW_T_R );   %D2��ÿ������ĸ��ʷֲ����

%ÿ�������µ��ʵĸ���
K_R = 100;                                     %��ʾǰ50������(�˴��ġ����ʡ�Ϊservice)
Sorted_P2_T_S = zeros( K_R , T_R );           %D2��ÿ��������service�ķֲ������D2_topic_serviceֵ�Ĵ�С(������)
Index_P2_T_S = zeros( K_R , T_R );            %D2��ÿ��������service���±�ֲ���D2_topic_service����service���±�(������)

for t=1:T_R
   [ temp1 , temp2 ] = sort( -W_T_R( : , t ) );
   Sorted_P2_T_S( : , t )  = ( full( -temp1( 1:K_R )) + BETA ) ./ ( repmat( sumW_T_R( t ) , K_R , 1 ));
   Index_P2_T_S( : , t )   = temp2( 1:K_R );
end

sNum_R = size(S_T_R,1);
P2_T_S = zeros( sNum_R , T_R );             %D2��ÿ��������service�ķֲ������D2_topic_serviceֵ�Ĵ�С��δ����
for t=1:T_R
    temp1 = W_T_R(:,t);
    P2_T_S( : , t )  = ( temp1 + BETA ) ./ ( repmat( sumW_T_R( t ) , sNum_R , 1 ));
end

P2_S_T = zeros(sNum_R, T_R);                %D2��ÿ���ĵ�(service)������ֲ������D2_service_topicֵ�ô�С
for i=1:sNum_R
    for t=1:T_R
        P2_S_T(i,t) = (S_T_R(i,t) + ALPHA)/(sum(S_T_R(i,:)) + ALPHA*T_R);
    end
end
%%
% ���������ļ�
save STEP3DATA Sorted_P2_T_S Index_P2_T_S P2_S_T P2_T_S WO_R W_T_R S_T_R Z_R K_R probtopic_R
%%

