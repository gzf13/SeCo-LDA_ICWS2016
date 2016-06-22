%%
% 12.23
% ����������Ŀ��ѡ��
% mashup - service

%%
% mashup-service
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
% �����Ե���������
T_N = 5:5:70;
LValue = zeros(1,length(T_N));      % ��Ȼ����ֵ
minD_JS = zeros(1,length(T_N));     % ƽ����С�������ƶ�
perplexity = zeros(1,length(T_N));  % �����


for t_i = 1:length(T_N)
   disp('**************************************************************')
   disp(['��ǰ������Ŀ:',num2str(T_N(t_i))]);
   disp('1. ��ʼLDAģ�ͼ���...');
   % LDAģ��
   [WS, DS] = importworddoccounts('step5_0_wordcount_msRelation.txt');    %����importworddoccounts�����������

    WO = service(ssum,1);
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
    T_R
    sumW_T_R = sum( W_T_R , 1 ) + BETA*W_R;
    probtopic_R = sumW_T_R / sum( sumW_T_R );   %D2��ÿ������ĸ��ʷֲ����

    %P_T_S
    sNum_R = size(S_T_R,1);
    P_T_S = zeros( W_R , T_R );             
    for t=1:T_R
        temp1 = W_T_R(:,t);
        P_T_S( : , t )  = ( temp1 + BETA ) ./ ( repmat( sumW_T_R( t ) , W_R , 1 ));
    end
    %P_M_T
    P_M_T = zeros(sNum_R, T_R);                
    for i=1:sNum_R
        P_M_T(i,:) = (S_T_R(i,:) + ALPHA)/(sum(S_T_R(i,:)) + ALPHA*T_R);
    end

    % ������Ȼ����ֵ
    disp('3. ��ʼ��Ȼ����ֵ����...')
    temLValue = 0;
    for i=1:msMNum_new
       % each doucment(mashup)
       word_index = find(msRelation(:,i));  %�ҳ�di���������е���
       for w_i = 1:length(word_index)
           n_wi = 1;
           temLValue = temLValue + n_wi * log(sum(P_T_S(word_index(w_i),:).*P_M_T(i,:)));     %sigma[z]P(w|z)P(z|d)
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
               D_JS(a,b) =  calJS(P_T_S(:,a),P_T_S(:,b));
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
    for i=1:msMNum_new
        % each document(mashup here)
        word_index = find(msRelation(:,i));  %�ҳ�di���������е���
        for w_i = 1:length(word_index)
           sigma_log_pwd = sigma_log_pwd + log(sum(P_T_S(word_index(w_i),:).*P_M_T(i,:)));     %sigma[z]P(w|z)P(z|d)
       end
    end
    sigma_Nd = sum(sum(msRelation));
    perplexity(t_i) = exp(-sigma_log_pwd/sigma_Nd);
end


save STEP5_0_mcitation minD_JS perplexity;
%%
% % ��������
% plot(T_N,LValue,'-b*');
% [~,index] = max(LValue);
% disp('******************************************************')
% disp(['���������Ϊ��',num2str(T_N(index))]);
% 
% 
% %%
% K_R = 10;                                     %��ʾǰ50������(�˴��ġ����ʡ�Ϊservice)
% Sorted_P2_T_S = zeros( K_R , T_R );           %D2��ÿ��������service�ķֲ������D2_topic_serviceֵ�Ĵ�С(������)
% Index_P2_T_S = zeros( K_R , T_R );            %D2��ÿ��������service���±�ֲ���D2_topic_service����service���±�(������)
% 
% for t=1:T_R
%    [ temp1 , temp2 ] = sort( -W_T_R( : , t ) );
%    Sorted_P2_T_S( : , t )  = ( full( -temp1( 1:K_R )) + BETA ) ./ ( repmat( sumW_T_R( t ) , K_R , 1 ));
%    Index_P2_T_S( : , t )   = temp2( 1:K_R );
% end


%%
% plot(T_N,perplexity,'--bs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','b');
% hold on
% plot(T_N,XX,'--rs','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','r');
% title('perplexity');
% legend('masup-service citation','service-service citation')