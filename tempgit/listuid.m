for i = 1:length(m)
    uuid_list(i).uuid=char(m(1).UUID.toString())
end      
%save('parralelisme_uuidlist.mat','uuid_list')