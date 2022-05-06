nvars = 4;
lb = [20,0,2,0.5];
ub = [70,0.2095,15,7];
initial_x = [ones(100,1)*60,ones(100,1)*0.11,ones(100,1)*10,ones(100,1)*6];

options = optimoptions(@particleswarm,'OutputFcn',@myout);

options.InitialSwarmMatrix = initial_x;
options.MaxIterations = 30;
options.SwarmSize = 5;
options.UseVectorized = true;
options.Display = 'iter';
%options = optimoptions(options,'PlotFcn',@pswplotbestf);

NN_pred = @(x)NN_p(x);

n = 0;
%eval_store = [];

while n < 2
inter_swarmfvals = [];
save(fullfile(pwd,'inter_swarmfvals'),'inter_swarmfvals');
[final_x,fval,exitflag,output] = particleswarm(NN_pred,nvars,lb,ub,options);

if fval>0
    n=n+1;
    load('inter_swarmfvals.mat');
    eval_store = [eval_store,inter_swarmfvals];    
end

end






function outputs = NN_p(inputs)
    NN = py.tensorflow.keras.models.load_model('evaluate_model');
    scaler_x = py.pickle.load(py.open('Xscaler.sav', 'rb'));
    scaler_y = py.pickle.load(py.open('Yscaler.sav', 'rb'));
    
    inputs_now = mat2np(inputs');
    inputs_now = inputs_now.reshape(int32(-1),int32(4));
    inputs_now = scaler_x.transform(inputs_now);
    
    outputs_now = NN.predict(inputs_now);
    outputs_now = scaler_y.inverse_transform(outputs_now);
    outputs_now = np2mat(outputs_now);
    outputs = outputs_now(:,2)./outputs_now(:,1);
    %outputs(:,3) = -outputs_now(2)+48;
end

function stop = myout(optimValues,state)
inter_swarmfvals = optimValues.swarmfvals;
saved_swarmfvals = load('inter_swarmfvals.mat');
inter_swarmfvals = [saved_swarmfvals.inter_swarmfvals;inter_swarmfvals];
stop = false; % This function does not stop the solver
save(fullfile(pwd,'inter_swarmfvals'),'inter_swarmfvals')
end