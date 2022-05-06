function outputs = ftir_NN_2D(inputs)
    NN = py.tensorflow.keras.models.load_model('evaluate_model');
    scaler_x = py.pickle.load(py.open('Xscaler.sav', 'rb'));
    scaler_y = py.pickle.load(py.open('Yscaler.sav', 'rb'));
    
    size_inputs = size(inputs);
    inputs_now = [repmat(50,size_inputs(1),1), inputs, repmat(7,size_inputs(1),1)];
    inputs_now = mat2np(inputs_now');
    inputs_now = inputs_now.reshape(int32(-1),int32(4));
    inputs_now = scaler_x.transform(inputs_now);
    
    outputs_now = NN.predict(inputs_now);
    outputs_now = scaler_y.inverse_transform(outputs_now);
    outputs_now = np2mat(outputs_now);
    outputs = outputs_now(:,2)./outputs_now(:,1);
    for i = 1:length(outputs)
       if outputs_now(i,1) <= 0
           a = 0.1;
       else
           a = outputs_now(i,1);
       end
       
       if outputs_now(i,2) <= 0
           b = 0.1;
       else
           b = outputs_now(i,2);
       end
       outputs(i) = b/a;
    end
    outputs = outputs.*(60*24*1e-3/(14*2.5*1e-6));
end


