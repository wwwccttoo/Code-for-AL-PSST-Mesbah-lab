function [share_samples,share_values,best_point,best_F,best_C] = ADMMBO(problem,opt)  
%% runs the ADMMBO algaorithm given initial points & budgets using EI

    fprintf('#######################################\n')
    fprintf('START OF "ADMMBO" ALGORTIHM\n');
    %% Preallocation
    ABSTOL=0.05; % ADMM Convergence In Practice Parameter
    RELTOL=0.05; %ADMM Convergence In Practice Parameter
    relaxp=0.5; % if >1 sets a faster convergence
    relaxd=1.;
    C_check=ones(length(problem.c),1);
    con=zeros(length(problem.c),1);
    mu=10; % Adapting rho Parameter
    thai=2; % Adapting rho Parameter
    Samples=[];
    tolC = 0.01;

    %% rebuildig the functions & evaluate the function values for an infeasible point
    if isfield(problem,'coeval')
        if problem.coeval == 1
            F = problem.coevalF;
            C = problem.coevalC;
            opt.f.coeval = 1;
            for i = 1:length(problem.c)
                opt.c{i}.coeval = 1;
            end
            share_samples = [];
            share_values = [];
        else
            F = problem.F;
            C = problem.C;
        end
    else
        problem.coeval = 0;
        F = problem.F;
        C = problem.C;
    end
    
    best_point = problem.InfPoint;
    F_prob = F(best_point);
    
    best_F = 1000 * F_prob(1) * sign(F_prob(1)); % a very positive value
    best_C = ones(length(problem.c),1);
    
    for s=1:length(problem.c)
        Cs_prob = C{s}(best_point);
        
        best_C(s) = Cs_prob(s+1)*1000 * sign(Cs_prob(s+1));
    end

    %% ADMM Main Loop
     for i=1:opt.ADMM.max_iters
        fprintf('#######################################\n')
        fprintf('Solving the Optimality subproblem at outer iteration %d of ADMMBO:\n',i);
        fprintf('    ******************************    \n')
        fprintf('Bayesian Optimizations inner iterations starts:\n')
        fprintf('    ******************************   \n')
        
        %% X-update 
        % building the Augmented Lagrangian function
        ybar=mean(cat(3,opt.f.y{:}),3);
        zbar=mean(cat(3,opt.f.z{:}),3);
        AL_X = @(X) F(X)+sum((length(problem.c)*opt.ADMM.rho/2)*(X-(zbar)+(ybar/opt.ADMM.rho)).^2,2);
        
        AL_part = @(X) sum((length(problem.c)*opt.ADMM.rho/2)*(X-(zbar)+(ybar/opt.ADMM.rho)).^2,2);

         % Optimizing the AUGMENTED LAGRANGIAN wrt X using BO
        if i~=1      
            if problem.coeval ~= 1
                opt.f.AL_evals = opt.f.true_evals+AL_part(opt.f.samples);
            else
                opt.f.AL_evals = share_values+AL_part(share_samples);
                opt.f.samples = share_samples;
            end
        end
        
        if i==1
            opt.f.max_iters=size(opt.f.initial_points,1)+(opt.f.step_iters);
        else
            opt.f.max_iters=size(opt.f.samples,1)+(opt.f.reduced_step_iters);
        end

        
                
        [xmin,~,T_u]=bayesopt(AL_X,opt.f); 
        XX(i,:)=xmin;
        
        true_evls = T_u.values - AL_part(T_u.samples);

        % Updating the Samples set based on X-step BO results
        oldSamplesNum=size(T_u.samples,1);
        if i==1    
            Samples=[Samples;T_u.samples];
        else
            Samples=[Samples;T_u.samples(oldSamplesNum+1:end,:)];
        end
        
        % updating x*_i for Z-update and Y-updates & gathering the
        % observed data
        opt.f.x=xmin;% Updating x* for Z-update and Y-updates
        %opt.f.initial_points=T_u.samples;% Updating initial X points/ this is CRAZY!
        opt.f.initial_points=[];
        opt.f.samples = T_u.samples;
        opt.f.true_evals = true_evls;
        if problem.coeval == 1
            share_values = opt.f.true_evals;
            share_samples = opt.f.samples;
        end
        %{
        disp(T_u.samples)
        disp(T_u.values)
        disp(opt.f.true_evals)
        disp(T_u.values(1)-delete_AL_part(T_u.samples(1,:)))
        disp(delete_AL_part(T_u.samples))
        %}
        
        if i==1
            iter_num = sprintf('%dst',i);
        else
            iter_num = sprintf('%dnd',i);
        end    
          
        % Updating the incumbent 
        if problem.coeval ~= 1
            [best_point,best_F,best_C] = incumbent_update(best_point,best_F,best_C,xmin,F,C,tolC,iter_num,'Optimality');
        else
            current_values_satisfy_cons = share_values(all(share_values(:,2:end)<=tolC,2),1);
            current_minvalues_satisfy_cons = min(current_values_satisfy_cons);
            [~, min_id_true] = ismember(current_minvalues_satisfy_cons, share_values(:,1), 'rows');
            xmin_undercon = share_samples(min_id_true,:);
            if isempty(xmin_undercon)
                [~, min_id] = ismember(xmin, share_samples, 'rows');
                [best_point,best_F,best_C] = incumbent_update(best_point,best_F,best_C,share_values,F,C,tolC,iter_num,'Optimality',min_id,xmin);
            else
                [best_point,best_F,best_C] = incumbent_update(best_point,best_F,best_C,share_values,F,C,tolC,iter_num,'Optimality',min_id_true,xmin_undercon);
            end
        end
        %% Z-update 
        fprintf('#######################################\n')
        fprintf('Solving the Feasibility subproblem at %s outer iteration of ADMMBO:\n',iter_num);
        fprintf('    ******************************    \n')
        fprintf('Bayesian Optimizations inner iterations starts:\n')
        fprintf('    ******************************   \n')
        
        % Checking if we have already satisfied the constraint's coonvergence criterion by C_check
        if C_check
            % Keeping track of old Z* to check (Z*(k+1)-Z*(k))
            zold=opt.f.z;
            
            for j=1:length(problem.c)
                
                % Adapting the max number of BO iterations according to ADMMBO's inner loop
                if i==1
                    opt.c{j}.max_iters=opt.c{j}.step_iters;
                else
                    opt.c{j}.max_iters=opt.c{j}.reduced_step_iters;  
                end
                
                if problem.coeval == 1
                    opt.c{j}.initial_points=[];
                    opt.c{j}.samples=share_samples;
                    opt.c{j}.true_evals=share_values;
                end
                % Optimizing the feasibility subproblem j^th  wrt Z using BO
                [zmin{j},min_id,T_h{j}] = bayesfeas(problem,opt,j); 
                
                if j ==1
                    subproblem = [ sprintf('%dst',j) ' Feasibility'];
                else
                    subproblem = [ sprintf('%dnd',j) ' Feasibility'];
                end
                
                
                if i==1    
                    Samples=[Samples;T_h{j}.samples];
                else
                    Samples=[Samples;T_h{j}.samples(oldSamplesNumC{j}+1:end,:)];
                end
                
                ZZ{j}(i,:)=opt.f.z{j};
                % Updating the Samples set based on Z-step BO results
                oldSamplesNumC{j}=size(T_h{j}.samples,1);
                opt.f.z{j}=zmin{j};% Updating z* for X-update and Y-updates
                %opt.c{j}.initial_points=T_h{j}.samples;% Updating initial Z points/ this is CRAZY!!!!!!
                if problem.coeval == 1
                    share_samples=T_h{j}.samples;
                    share_values=T_h{j}.values;
                else
                    opt.c{j}.initial_points=[];
                    opt.c{j}.samples=T_h{j}.samples;
                    opt.c{j}.true_evals=T_h{j}.values;
                end
                % Updating the incumbent 
                if problem.coeval ~= 1
                    [best_point,best_F,best_C] = incumbent_update(best_point,best_F,best_C,zmin{j},F,C,tolC,iter_num,subproblem);
                else
                    current_values_satisfy_cons = share_values(all(share_values(:,2:end)<=tolC,2),1);
                    current_minvalues_satisfy_cons = min(current_values_satisfy_cons);
                    [~, min_id_true] = ismember(current_minvalues_satisfy_cons, share_values(:,1), 'rows');
                    xmin_undercon = share_samples(min_id_true,:);
                    if isempty(xmin_undercon)
                        [~, min_id] = ismember(zmin{j}, share_samples, 'rows');
                        [best_point,best_F,best_C] = incumbent_update(best_point,best_F,best_C,share_values,F,C,tolC,iter_num,'Optimality',min_id,zmin{j});
                    else
                        [best_point,best_F,best_C] = incumbent_update(best_point,best_F,best_C,share_values,F,C,tolC,iter_num,'Optimality',min_id_true,xmin_undercon);
                    end
                end
                %% Y-update 
                ymin = opt.f.y{j} +opt.ADMM.rho*(xmin - zmin{j});
                opt.f.y{j}=ymin;
                clear  ymin
            end
        end
        %% Check the termination Condition
        for j=1:length(problem.c) 
             history.r_norm{j}= norm(opt.f.x-opt.f.z{j});
             history.s_norm{j}=norm(-opt.ADMM.rho*(opt.f.z{j}- zold{j}));
             history.eps_pri{j} = 0.1*sqrt(opt.f.dims)*ABSTOL +relaxp*RELTOL*max(norm(opt.f.x), norm(-opt.f.z{j}));
             history.eps_dual{j}= 0.3*sqrt(opt.f.dims)*ABSTOL + relaxd*RELTOL*norm(opt.ADMM.rho*opt.f.y{j});
             history.r_norm{j}
             history.s_norm{j}
             history.eps_pri{j}
             history.eps_dual{j}
             if (history.r_norm{j} < history.eps_pri{j} && ...
              history.s_norm{j} < history.eps_dual{j})
              con(j)=1;
             end
        end
         
      % checking the ADMM convergence
     if con
       fprintf("It takes %d ADMM iterations To converge.\n ",i);
         %break;
     end
     
     % updating the penalty parameter
     if mean(cell2mat(history.r_norm)) > mu*mean(cell2mat(history.s_norm))
        opt.ADMM.rho= opt.ADMM.rho*thai;
        fprintf("rho is increased\n");
     elseif  mu* mean(cell2mat(history.r_norm)) < mean(cell2mat(history.s_norm))
        opt.ADMM.rho= opt.ADMM.rho/thai;
        fprintf("rho is decreased\n");
     else
        opt.ADMM.rho= opt.ADMM.rho;
        fprintf("rho is unchanged\n");
     end
     
        
       
     end
end

function[best_point_updated,best_F_updated,best_C_updated] = incumbent_update(best_point,best_F,best_C,current_point,F,C,tolC,iter_num,subproblem,index,current_min)
    % this function updates and reports the incumbent and its objective value
    if nargin < 10
        current_F = F(current_point);
        current_F = current_F(1);
        current_C = ones(length(C),1);
        for s=1:length(C)
            inter_C = C{s}(current_point);
            current_C(s) = inter_C(s+1);
        end
    else
        current_F = current_point(index,1); % this current_point will be evaluations
        for s=1:length(C)
            current_C(s) = current_point(index,s+1);
        end
        current_point = current_min;
    end
    
    if sum(current_C > tolC) && sum(best_C >tolC)
        best_point_updated = best_point;
        best_F_updated = best_F; 
        best_C_updated = best_C;
        fprintf('    ******************************    \n')
        fprintf('No feasible point is found yet after %s ADMMBO iteration during %s subproblem! \n',iter_num,subproblem);
        fprintf('    ******************************   \n')
        
    elseif  sum(current_C>tolC) == 0 && current_F <= best_F
       best_point_updated = current_point;
       best_F_updated = current_F;
       best_C_updated = current_C;
       fprintf('    ******************************    \n')
       fprintf('The incumbent is updated & the bestfeasible oberved value after %s ADMMBO iteration during %s subproblem is %f.\n',iter_num,subproblem,best_F_updated);
       fprintf('    ******************************   \n')
    else
       best_point_updated = best_point;
       best_F_updated = best_F; 
       best_C_updated = best_C;
       fprintf('    ******************************    \n')
       fprintf('The incumbent is NOT updated at %s ADMMBO iteration during %s subproblem. \n',iter_num,subproblem);
       fprintf('    ******************************   \n')
       
         
    end
    
end
