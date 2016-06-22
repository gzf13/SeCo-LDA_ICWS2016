%%
i = 35;
i1 = leftServiceSet(Index_P2_T_S(1,i));
i2 = leftServiceSet(Index_P2_T_S(2,i));
i3 = leftServiceSet(Index_P2_T_S(3,i)); 
sum=0;
for j=1:mNum
    records = useRecord(:,j);
    if (records(i1)~=0 && records(i2)~=0)
        disp(char(mashup(j,1)));
        sum=sum+1;
    end
end