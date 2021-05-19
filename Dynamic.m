%% Change subName here to match the names of the edf files
subName = {''};
N=32;  % Noodes, which equals to the Channel number
T=4;   % Slices, each trial is 4s in length
% Each multi-slice network (of each trial) has 4 subnetworks each of
% which based on one-second signals
fs=128;
%% Initialize Variables to be further used
% a: alpha band, b: beta band, g: gamma band
% r: resting state, rand: random network (null model), otherwise task state network
Qa=[]; Qb = []; Qg = [];Qa_rand = []; Qb_rand = []; Qg_rand = [];
Qa_r = []; Qb_r = []; Qg_r = []; Qa_rand2 = []; Qb_rand2 = []; Qg_rand2 = [];
%%
for subInd = 1 : length(subName)
    %% Extract Blocks based on events
    signal = cell2mat(struct2cell(load(append(subName{subInd},'.mat'))));
    eventInf = xlsread(append(subName{subInd},'.csv'));
    eventTyp = eventInf(:,4);
    eventTime = round(fs.*eventInf(:,1)); % Convert Time to index
    
    % Search based on event types (may need to change values)
    % Search for statring points of the blocks
    eventLoc = round(find(eventTyp==1));
    eventLoc = [eventLoc; find(eventTyp==2)];
    rstLoc = round(find(eventTyp==3));
    % Search for endpoints of the blocks
    task_all = {};
    for i = 1 : length(eventLoc)-1
        task_all{i} = signal(:,eventTime(eventLoc(i)):eventTime(eventLoc(i+1)));
    end
    task_all(length(eventLoc)) = signal(eventLoc(end):end);
    for blockInd = 1 : length(Te)
        task_full = task_all{blockInd};
        task_full = task_full(channInclud,:);
        [chann_size,blockSize] = size(task_full);
        tPoint = 1:fs:blockSize;
        %% Task Trials
        for timeInd = 1 : 8
            task = task_full(:,tPoint(timeInd):tPoint(timeInd+1));
            cohMap_task_alpha = zeros(chann_size,chann_size);
            cohMap_task_beta = zeros(chann_size,chann_size);
            cohMap_task_gamma = zeros(chann_size,chann_size);
            [~,~,f] = bst_cohn(task(1,:), task(2,:), fs, 3, 0.5, 'icohere2019');
            % Find indexes of the frequency bands 8-13 Hz
            [~, ind1_alpha] = min(abs(f-8));
            [~, ind2_alpha] = min(abs(f-13));
            % Find indexes of the frequency bands 14-30 Hz
            [~, ind1_beta] = min(abs(f-14));
            [~, ind2_beta] = min(abs(f-30));
            % Find indexes of the frequency bands 30-40 Hz
            [~, ind1] = min(abs(f-30));
            [~, ind2] = min(abs(f-40));
            % Calculated coherences
            for i = 1 : chann_size
                for j = 1 : chann_size
                    coh = bst_cohn(task(i,:), task(j,:), fs, 3, 0.5, 'icohere2019',0, [], 100);
                    coh = squeeze(coh);
                    cohMap_task_alpha(i,j) = mean(coh(ind1_alpha:ind2_alpha));
                    cohMap_task_beta(i,j) = mean(coh(ind1_beta:ind2_beta));
                    cohMap_task_gamma(i,j) = mean(coh(ind1:ind2));
                end
            end
            
            if timeInd <= 4
                A_r{timeInd} = weight_conversion(ones(36,36)-cohMap_task_alpha,'autofix');
                A_rand2{timeInd} = null_model_und_sign(A_r{timeInd});
                B_r{timeInd} = weight_conversion(ones(36,36)-cohMap_task_beta,'autofix');
                B_rand2{timeInd} = null_model_und_sign(B_r{timeInd});
                G_r{timeInd} = weight_conversion(ones(36,36)-cohMap_task_gamma,'autofix');
                G_rand2{timeInd} = null_model_und_sign(G_r{timeInd});
            end
            if timeInd >4
                A{timeInd-4} = weight_conversion(ones(36,36)-cohMap_task_alpha,'autofix');
                A_rand{timeInd-4} = null_model_und_sign(A{timeInd-4});
                B{timeInd-4} = weight_conversion(ones(36,36)-cohMap_task_beta,'autofix');
                B_rand{timeInd-4} = null_model_und_sign(B{timeInd-4});
                G{timeInd-4} = weight_conversion(ones(36,36)-cohMap_task_gamma,'autofix');
                G_rand{timeInd-4} = null_model_und_sign(G{timeInd-4});
            end
        end
        
        %% Dynamic Modularity
        B1=spalloc(N*T,N*T,N*N*T+2*N*T);
        B2=spalloc(N*T,N*T,N*N*T+2*N*T);
        B3=spalloc(N*T,N*T,N*N*T+2*N*T);
        B1_r=spalloc(N*T,N*T,N*N*T+2*N*T);
        B2_r=spalloc(N*T,N*T,N*N*T+2*N*T);
        B3_r=spalloc(N*T,N*T,N*N*T+2*N*T);
        B1_rand=spalloc(N*T,N*T,N*N*T+2*N*T);
        B2_rand=spalloc(N*T,N*T,N*N*T+2*N*T);
        B3_rand=spalloc(N*T,N*T,N*N*T+2*N*T);
        B1_rand2=spalloc(N*T,N*T,N*N*T+2*N*T);
        B2_rand2=spalloc(N*T,N*T,N*N*T+2*N*T);
        B3_rand2=spalloc(N*T,N*T,N*N*T+2*N*T);
        twomu1=0;twomu2=0;twomu3=0;
        gamma = 1;omega = 1;
        for s=1:T
            k1=sum(A{s});k1r=sum(A_r{s});
            k2=sum(B{s});k2r=sum(B_r{s});
            k3=sum(G{s});k3r=sum(G_r{s});
            k1_rand=sum(A_rand{s});k1_rand2=sum(A_rand2{s});
            k2_rand=sum(B_rand{s});k2_rand2=sum(B_rand2{s});
            k3_rand=sum(G_rand{s});k3_rand2=sum(G_rand2{s});
            twom1=sum(k1);twom2=sum(k2);twom3=sum(k3);
            twomu1=twomu1+twom1;twomu2=twomu2+twom2;twomu3=twomu3+twom3;
            indx=[1:N]+(s-1)*N;
            B1(indx,indx)=A{s}-gamma*k1'*k1/twom1;
            B1_r(indx,indx)=A_r{s}-gamma*k1r'*k1r/twom1;
            B1_rand(indx,indx)=A_rand{s}-gamma*k1_rand'*k1_rand/twom1;
            B1_rand2(indx,indx)=A_rand2{s}-gamma*k1_rand2'*k1_rand2/twom1;
            B2(indx,indx)=B{s}-gamma*k2'*k2/twom2;
            B2_r(indx,indx)=B_r{s}-gamma*k2r'*k2r/twom2;
            B2_rand(indx,indx)=B_rand{s}-gamma*k2_rand'*k2_rand/twom2;
            B2_rand2(indx,indx)=B_rand2{s}-gamma*k2_rand2'*k2_rand2/twom2;
            B3(indx,indx)=G{s}-gamma*k3'*k3/twom3;
            B3_r(indx,indx)=G_r{s}-gamma*k3r'*k3r/twom3;
            B3_rand(indx,indx)=G_rand{s}-gamma*k3_rand'*k3_rand/twom3;
            B3_rand2(indx,indx)=G_rand2{s}-gamma*k3_rand2'*k3_rand2/twom3;
        end
        twomu1=twomu1+2*omega*N*(T-1);twomu2=twomu2+2*omega*N*(T-1);twomu3=twomu3+2*omega*N*(T-1);
        B1 = B1 + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B1_r = B1_r + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B1_rand = B1_rand + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B1_rand2 = B1_rand2 + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B2 = B2 + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B2_r = B2_r + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B2_rand = B2_rand + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B2_rand2 = B2_rand2 + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B3 = B3 + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B3_r = B3_r + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B3_rand = B3_rand + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        B3_rand2 = B3_rand2 + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        [S1,Q1] = genlouvain(B1);
        [~,Q1_r] = genlouvain(B1_r);
        [~,Q1_rand] = genlouvain(B1_rand);
        [~,Q1_rand2] = genlouvain(B1_rand2);
        [S2,Q2] = genlouvain(B2);
        [~,Q2_r] = genlouvain(B2_r);
        [~,Q2_rand] = genlouvain(B2_rand);
        [~,Q2_rand2] = genlouvain(B2_rand2);
        [S3,Q3] = genlouvain(B3);
        [~,Q3_r] = genlouvain(B3_r);
        [~,Q3_rand] = genlouvain(B3_rand);
        [~,Q3_rand2] = genlouvain(B3_rand2);
        Qa = [Qa Q1/twomu1];Qb = [Qb Q2/twomu2];Qg = [Qg Q3/twomu3];
        Qa_r = [Qa_r Q1_r/twomu1];Qb_r = [Qb_r Q2_r/twomu2];Qg_r = [Qg_r Q3_r/twomu3];
        Qa_rand = [Qa_rand Q1_rand/twomu1];Qb_rand = [Qb_rand Q2_rand/twomu2];Qg_rand = [Qg_rand Q3_rand/twomu3];
        Qa_rand2 = [Qa_rand2 Q1_rand2/twomu1];Qb_rand2 = [Qb_rand2 Q2_rand2/twomu2];Qg_rand2 = [Qg_rand2 Q3_rand2/twomu3];
    end
    % Report Progress
    subInd
end