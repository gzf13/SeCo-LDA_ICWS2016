%%
% 12.23
% ����������Ŀ��ѡ��

%%
% service-service


load('../Record.mat');
load('STEP2DATA.mat');

%%
% �����Ե���������
T_N = 5:5:70;
LValue = zeros(1,length(T_N));
minD_JS = zeros(1,length(T_N));     % ƽ����С�������ƶ�
perplexity = zeros(1,length(T_N));  % �����


for t_i = 1:length(T_N)
   disp('**************************************************************')
   disp(['��ǰ������Ŀ:',num2str(T_N(t_i))]);
   disp('1. ��ʼLDAģ�ͼ���...');
   % LDAģ��
   [WS, DS] = importworddoccounts('step2_wordcount_ssRelationLeft.txt');    %����importworddoccounts�����������

    WO = service(leftServiceSet,1);
    WO_R = WO';

    % Set the number of topics
    % ��t_i��ʵ��
    T = T_N(t_i);

    % Set the hyperparameters
    BETA=0.01;
    ALPHA=1;

    % The number of iterations
    N = 400; 
    % The random seed
    SEED = 3;
    % What output to show (0=no output; 1=iterations; 2=all output)
    OUTPUT = 0;
    % This function might need a few minutes to finish
    tic
    %[WP,DP,Z]
    [ W_T_R,S_T_R,Z_R ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
    toc
    
    % ����P_S_T && P_T_S
    disp('2. ��ʼ��ظ��ʷֲ�����...');
    
    %ÿ������ĸ���
    W_R = size( W_T_R , 1 );                    %���ʵ�������D2��Ϊservice������
    T_R = size( W_T_R , 2 );                    %���������
    %T_R
    sumW_T_R = sum( W_T_R , 1 ) + BETA*W_R;
    probtopic_R = sumW_T_R / sum( sumW_T_R );   %D2��ÿ������ĸ��ʷֲ����

    %P_T_S
    sNum_R = size(S_T_R,1);
    P2_T_S = zeros( sNum_R , T_R );             %D2��ÿ��������service�ķֲ������D2_topic_serviceֵ�Ĵ�С��δ����
    for t=1:T_R
        temp1 = W_T_R(:,t);
        P2_T_S( : , t )  = ( temp1 + BETA ) ./ ( repmat( sumW_T_R( t ) , sNum_R , 1 ));
    end
    %P_S_T
    P2_S_T = zeros(sNum_R, T_R);                %D2��ÿ���ĵ�(service)������ֲ������D2_service_topicֵ�ô�С
    for i=1:sNum_R
        P2_S_T(i,:) = (S_T_R(i,:) + ALPHA)/(sum(S_T_R(i,:)) + ALPHA*T_R);
    end

    % ������Ȼ����ֵ
    disp('3. ��ʼ��Ȼ����ֵ����...')
    temLValue = 0;
    for i=1:leftServiceNum
       % each doucment(service here)
       word_index = find(ssRelationLeft(i,:));  %�ҳ�di���������е���
       for w_i = 1:length(word_index)
           n_wi = ssRelationLeft(i,word_index(w_i));
           temLValue = temLValue + n_wi * log(sum(P2_T_S(word_index(w_i),:).*P2_S_T(i,:)));     %sigma[z]P(w|z)P(z|d)
       end
    end
    LValue(t_i) = temLValue;
    disp(['��Ȼ����ֵΪ:',num2str(temLValue)]);
    
    % ����ƽ����С�������ƶ�
    disp('4. ��ʼƽ����С�������ƶȼ���...');
    D_JS = zeros(T,T);
    temAverJS = 0;
    for a = 1:T
        for b = 1: T
            if (a==b) 
                D_JS(a,b) = 9999;
            else
               D_JS(a,b) =  calJS(P2_T_S(:,a),P2_T_S(:,b));
            end
        end
    end
    for i = 1:T
        temAverJS = temAverJS + min(D_JS(i,:));
    end
    temAverJS = temAverJS/T;
    minD_JS(t_i) = temAverJS;
    
    % ���������
    disp('5. ��ʼ����ȼ���...');
    sigma_log_pwd = 0;
    for i=1:leftServiceNum
        % each document(service here)
        word_index = find(ssRelationLeft(i,:));  %�ҳ�di���������е���
        for w_i = 1:length(word_index)
           sigma_log_pwd = sigma_log_pwd + log(sum(P2_T_S(word_index(w_i),:).*P2_S_T(i,:)));     %sigma[z]P(w|z)P(z|d)
       end
    end
    sigma_Nd = sum(sum(ssRelationLeft));
    perplexity(t_i) = exp(-sigma_log_pwd/sigma_Nd);
end

save STEP5_0_scitation minD_JS perplexity;

% % ��������
% plot(T_N,LValue,'-b*');
% [~,index] = max(LValue);
% disp('******************************************************')
% disp(['���������Ϊ��',num2str(T_N(index))]);