% subName = {''};
p_beta_right = 0; p_beta_left = 0;
p_alpha_right = 0; p_alpha_left = 0;
laterality1 = zeros(20,1);laterality2 = zeros(20,1);
Fs = Fs;
%%
for subInd = 1 : 20
    %% Extract Blocks
    signal = cell2mat(struct2cell(load(append(subName{subInd},'.mat'))));
    [~,header] = sload(append(subName{subInd},'.edf'));
    eventTyp = header.EVENT.TYP;
    % Change this part to search for left/right hand blocks
%     eventLoc = find(eventTyp==3); 
%     eventLoc = find(eventTyp==2);
   
    rstLoc = find(eventTyp==1);
    Ts = header.EVENT.POS(eventLoc);
    Te = header.EVENT.POS(rstLoc);
    task_all = {};
    % change 9 and 13 to the corresponding number of C3 and C4
    for i = 1 : length(Ts)
        task_all_l{i} = signal(9,Ts(i)+0.5*Fs:Ts(i) + 3.5*Fs);
        task_all_r{i} = signal(13,Ts(i)+0.5*Fs:Ts(i) + 3.5*Fs);
    end
    rst_all ={};
    Te = Ts;
    Ts = header.EVENT.POS(eventLoc-1);
    for i = 1 : length(Te)
        rst_all_l{i} = signal(9,Ts(i)+0.5*Fs:Ts(i)+3.5*Fs);
        rst_all_r{i} = signal(13,Ts(i)+0.5*Fs:Ts(i)+3.5*Fs);
    end 
    p_ERD_right = getERD(task_all_r,rst_all_r);
    p_ERD_left = getERD(task_all_l,rst_all_r);
    p_ERD_right = p_ERD_right-min([p_ERD_right,p_ERD_left]);
    p_ERD_left = p_ERD_left-min([p_ERD_right,p_ERD_left]);
    % Left hand
%          laterality1(subInd) = (mean(p_ERD_right)-mean(p_ERD_left))/(mean(p_ERD_right)+mean(p_ERD_left));
    % Right hand
   laterality2(subInd) = (mean(p_ERD_left)-mean(p_ERD_right))/(mean(p_ERD_right)+mean(p_ERD_left));
end