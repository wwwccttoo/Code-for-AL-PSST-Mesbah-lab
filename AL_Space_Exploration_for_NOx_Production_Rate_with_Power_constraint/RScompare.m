RS_final_list = [];
BO_final_list = [];
load('50ADMMBOruns');
RS_param_store = [];
for i = 1:50
   RS_thistime = [55, 0.05, 12.0, 5.0;rand(271,1)*50+20,rand(271,1)*0.2095,rand(271,1)*5+2,rand(271,1)*6.5+0.5];
   RS_param_store = [RS_param_store;RS_thistime];
end

RS_value_store = ftir_NN(RS_param_store,0);

for i = 1:50
    RS_thislist = [];
    BO_thislist = [];
    current_total = RS_value_store((i-1)*272+1:i*272,:);
    current_RS_best = current_total(1);
    current_BO_best = Values{i}(1,1);
    for j = 1:271
        RS_NOx = current_total(j,:);
        if current_RS_best >= RS_NOx(1) && RS_NOx(2) <=50
            current_RS_best = RS_NOx(1);
        end
        RS_thislist = [RS_thislist;current_RS_best];
        
        
        
        if current_BO_best >= Values{i}(j,1) && Values{i}(j,2)<=0
            current_BO_best = Values{i}(j,1);
        end
        BO_thislist = [BO_thislist;current_BO_best];
    end
    RS_final_list = [RS_final_list,RS_thislist];
    BO_final_list = [BO_final_list,BO_thislist];
end


