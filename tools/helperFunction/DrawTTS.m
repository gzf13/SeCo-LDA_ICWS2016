function DrawTTS(TTS,average_TTS,i,totalT,sigma)
%%
% ���Ƶ�i�������ʱ����������

%%
temT = TTS(i,:);
[~, x0] = find(temT);
y0 = temT(x0);

%%
% ����ԭʼ���ݵ�

% 1. ֱ�ӻ��
figure();
plot(x0,y0,'rs');
hold on;
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
title(i);
rTTS = rTTS/max(rTTS)*max(temT);
plot(rTTS,'b');                                             %���Ʋ�ֵƽ�����ʱ����������
plot([average_TTS(i) average_TTS(i)],[0,max(temT)*2],'g');  %���ƾ�ֵ
axis([0 3800 0 max(temT)*1.2]);                             %���������᷶Χ

end
