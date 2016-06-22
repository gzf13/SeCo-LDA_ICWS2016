%%
%experiment3
%edit by gzf 1219
%STEP3:Relation-based LDA

%%
% Content-based LDA
% 参数准备

load('../Record.mat');
load('STEP2DATA.mat');
[WS, DS] = importworddoccounts('step2_wordcount_ssRelationLeft.txt');    %利用importworddoccounts读入所需参数

%使用清洗过的数据

%WO为service列表
%选取存在引用关系的service作为词典
%这些service在总的service下的编号记录在leftServiceSet中，数量为leftServiceNum

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
%预备工作二：整理D2相关数据
%DIC: WO_R
%[WP,DP,Z]
%[ W_T_R,S_T_R,Z_R ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
%D2:service-topic-service

%每个主题的概率
W_R = size( W_T_R , 1 );                    %单词的数量，D2中为service的数量
T_R = size( W_T_R , 2 );                    %主题的数量，50
sumW_T_R = sum( W_T_R , 1 ) + BETA*W_R;
probtopic_R = sumW_T_R / sum( sumW_T_R );   %D2，每个主题的概率分布情况

%每个主题下单词的概率
K_R = 100;                                     %显示前50个单词(此处的“单词”为service)
Sorted_P2_T_S = zeros( K_R , T_R );           %D2，每个主题中service的分布情况，D2_topic_service值的大小(已排序)
Index_P2_T_S = zeros( K_R , T_R );            %D2，每个主题中service的下标分布，D2_topic_service具体service的下标(已排序)

for t=1:T_R
   [ temp1 , temp2 ] = sort( -W_T_R( : , t ) );
   Sorted_P2_T_S( : , t )  = ( full( -temp1( 1:K_R )) + BETA ) ./ ( repmat( sumW_T_R( t ) , K_R , 1 ));
   Index_P2_T_S( : , t )   = temp2( 1:K_R );
end

sNum_R = size(S_T_R,1);
P2_T_S = zeros( sNum_R , T_R );             %D2，每个主题中service的分布情况，D2_topic_service值的大小（未排序）
for t=1:T_R
    temp1 = W_T_R(:,t);
    P2_T_S( : , t )  = ( temp1 + BETA ) ./ ( repmat( sumW_T_R( t ) , sNum_R , 1 ));
end

P2_S_T = zeros(sNum_R, T_R);                %D2，每个文档(service)的主题分布情况，D2_service_topic值得大小
for i=1:sNum_R
    for t=1:T_R
        P2_S_T(i,t) = (S_T_R(i,t) + ALPHA)/(sum(S_T_R(i,:)) + ALPHA*T_R);
    end
end
%%
% 保存结果到文件
save STEP3DATA Sorted_P2_T_S Index_P2_T_S P2_S_T P2_T_S WO_R W_T_R S_T_R Z_R K_R probtopic_R
%%

