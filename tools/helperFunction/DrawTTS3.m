function DrawTTS3(TTS,average_TTS,i,i2,totalT,sigma,Index_P2_T_S,Sorted_P2_T_S,service_date,minDate)
%%
% ���Ƶ�i�������ʱ����������



%%
%----------------------------------------------------------------------
% figure();
% for j=1:100
%     plot(service_date(Index_P2_T_S(j,i)),Sorted_P2_T_S(j,i),'ks','MarkerFaceColor','b');
%     hold on;
% end
%----------------------------------------------------------------------
%%
temT = TTS(i,:);
[~, x0] = find(temT);
y0 = temT(x0);

%%
% ����ԭʼ���ݵ�

% 1. ֱ�ӻ��

% 2. ����Ȩ�ش�С���ò�ͬ��makersize
% red,cubic
% for j=1:length(x0)
%     plot(x0(j),y0(j),'rs','markersize',20*y0(j)/max(y0));
%     hold on;
% end

% 3. ��ֵ����
rTTS = zeros(1,totalT);
for j=1:totalT
    tem = 0;
    for k=1:totalT
        % norm N(0,sigma)
        tem = tem + normpdf(k-j,0,sigma)/normpdf(0,0,sigma)*temT(k);
    end
    rTTS(j) = tem;
    j
end
rTTS = rTTS/sum(rTTS)*120;
plot(rTTS,'b');                                             %���Ʋ�ֵƽ�����ʱ����������
plot([average_TTS(i) average_TTS(i)],[0,max(temT)*2],'b');  %���ƾ�ֵ
axis([0 3800 0 max(temT)*1.2]);                             %���������᷶Χ

%%
% ���Ƶڶ�������
temT = TTS(i2,:);
[~, x02] = find(temT);
y02 = temT(x02);
%%
% ����ԭʼ���ݵ�

% 1. ֱ�ӻ��

% 2. ����Ȩ�ش�С���ò�ͬ��makersize
% red,cubic
% for j=1:length(x0)
%     plot(x0(j),y0(j),'rs','markersize',20*y0(j)/max(y0));
%     hold on;
% end

% 3. ��ֵ����
rTTS2 = zeros(1,totalT);
for j=1:totalT
    tem = 0;
    for k=1:totalT
        % norm N(0,sigma)
        tem = tem + normpdf(k-j,0,sigma)/normpdf(0,0,sigma)*temT(k);
    end
    rTTS2(j) = tem;
    j
end

rTTS2 = rTTS2/sum(rTTS2)*120;

%% ��ͼ
figure();
for j=1:100
    serviceIndex = Index_P2_T_S(j,i);
    serviceTime = datenum(service_date(serviceIndex,1),service_date(serviceIndex,2),service_date(serviceIndex,3)) - datenum(minDate(1),minDate(2),minDate(3)) + 1;
    plot(serviceTime,Sorted_P2_T_S(j,i),'ks','MarkerFaceColor','b');
    hold on;
end
%plot(x0,y0,'ks','MarkerFaceColor','b');
hold on;

f1=plot(rTTS*1.15,'--b','LineWidth',1.4);                                             %���Ʋ�ֵƽ�����ʱ����������
plot([average_TTS(i) average_TTS(i)],[0,max(temT)*3],'b');  %���ƾ�ֵ
axis([0 3800 0 max(temT)*1.3]);                             %���������᷶Χ

%plot(x02,y02,'ks','MarkerFaceColor','r');
for j=1:100
    serviceIndex = Index_P2_T_S(j,i2);
    serviceTime = datenum(service_date(serviceIndex,1),service_date(serviceIndex,2),service_date(serviceIndex,3)) - datenum(minDate(1),minDate(2),minDate(3)) + 1;
    plot(serviceTime,Sorted_P2_T_S(j,i2),'ks','MarkerFaceColor','r');
end

f2=plot(rTTS2,'--r','LineWidth',1.4);                                             %���Ʋ�ֵƽ�����ʱ����������
plot([average_TTS(i2) average_TTS(i2)],[0,max(temT)*3],'r');  %���ƾ�ֵ

%���ú�������ʾ����
set(gca,'XTick',120:365:3600);
set(gca,'XTickLabel',{'2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'});

str1 = strcat('Topic',num2str(i));
str2 = strcat('Topic',num2str(i2));
legend([f1,f2],str1,str2);

xlabel('Time');
ylabel('P(z|t)');

end
