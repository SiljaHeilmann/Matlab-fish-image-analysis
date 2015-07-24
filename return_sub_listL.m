function  L  = return_sub_listL(L,listOfIndex)

L_final = zeros(size(L));

for i=1:length(listOfIndex)
    index = listOfIndex(i);
    
    L_index = L;
    L_index(L>index)=0;
    L_index(L<index)=0;
    
    L_final = L_final + L_index;
    
end

L = L_final;

end

