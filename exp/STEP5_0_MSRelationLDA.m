%%
%experiment3
%edit by gzf 1220
%STEP5_0 MSRelation LDA (ֱ����m-s���ù�ϵ��LDA)

%%
% Content-based LDA
% ����׼��

load('../Record.mat');

%%
% ��ȥ����service����Ϊ0��1��mashup
msum = find(sum(useRecord)>1);
msMNum_new = length(msum);
msRelation = useRecord(:,msum);

% �ٳ�ȥû����mashup�г��ֹ���service
ssum = find(sum(msRelation,2)>0);
msSNum_new = length(ssum);
msRelation = msRelation(ssum,:);

%%
% ����ļ��Թ�LDA����ʹ��
fid = fopen('step5_0_wordcount_msRelation.txt','wt');
for i=1:1:msMNum_new
    %Mashup i
    temIndex = find(msRelation(:,i));
    temL = length(temIndex);
    for j=1:1:temL
        %������ļ�
        fprintf(fid,'%d %d %d\n',i,temIndex(j),1);
    end
    i
end
fclose(fid);

%%

[WS, DS] = importworddoccounts('step5_0_wordcount_msRelation.txt');    %����importworddoccounts�����������

%ʹ����ϴ��������

%WOΪservice�б�
%ѡȡ�������ù�ϵ��service��Ϊ�ʵ�

WO = service(ssum,1);
WO_MSR = WO';

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
[ W_T_MSR,S_T_MSR,Z_MSR ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
toc

%%
%Ԥ��������������D2�������
%DIC: WO_R
%[WP,DP,Z]
%[ W_T_R,S_T_R,Z_R ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
%D2:service-topic-service

%ÿ������ĸ���
W_MSR = size( W_T_MSR , 1 );                    %���ʵ�������D2��Ϊservice������
T_MSR = size( W_T_MSR , 2 );                    %�����������50
sumW_T_MSR = sum( W_T_MSR , 1 ) + BETA*W_MSR;
probtopic_MSR = sumW_T_MSR / sum( sumW_T_MSR );   %D2��ÿ������ĸ��ʷֲ����

%ÿ�������µ��ʵĸ���
K_R = 100;                                     %��ʾǰ50������(�˴��ġ����ʡ�Ϊservice)
Sorted_P3_T_S = zeros( K_R , T_MSR );           %D2��ÿ��������service�ķֲ������D2_topic_serviceֵ�Ĵ�С(������)
Index_P3_T_S = zeros( K_R , T_MSR );            %D2��ÿ��������service���±�ֲ���D2_topic_service����service���±�(������)

for t=1:T_MSR
   [ temp1 , temp2 ] = sort( -W_T_MSR( : , t ) );
   Sorted_P3_T_S( : , t )  = ( full( -temp1( 1:K_R )) + BETA ) ./ ( repmat( sumW_T_MSR( t ) , K_R , 1 ));
   Index_P3_T_S( : , t )   = temp2( 1:K_R );
end

sNum_R = size(W_T_MSR,1);
P3_T_S = zeros( sNum_R , T_MSR );             %D2��ÿ��������service�ķֲ������D3_topic_serviceֵ�Ĵ�С��δ����
for t=1:T_MSR
    temp1 = W_T_MSR(:,t);
    P3_T_S( : , t )  = ( temp1 + BETA ) / sumW_T_MSR( t );
end

P3_M_T = zeros(size(S_T_MSR,1), T_MSR);                %D2��ÿ���ĵ�(service)������ֲ������D3_mashup_topicֵ�ô�С
for i=1:size(S_T_MSR,1)
    for t=1:T_MSR
        P3_M_T(i,t) = (S_T_MSR(i,t) + ALPHA)/(sum(S_T_MSR(i,:)) + ALPHA*T_MSR);
    end
end
%%
% ���������ļ�
save STEP5_0_DATA Sorted_P3_T_S Index_P3_T_S P3_M_T P3_T_S WO_MSR W_T_MSR S_T_MSR Z_MSR K_R probtopic_MSR msum ssum msMNum_new msSNum_new msRelation T_MSR
%%

