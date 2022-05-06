function outputs = ftir_NN(inputs,param)
    NN = py.tensorflow.keras.models.load_model('evaluate_model');
    scaler_x = py.pickle.load(py.open('Xscaler.sav', 'rb'));
    scaler_y = py.pickle.load(py.open('Yscaler.sav', 'rb'));
    
    inputs_now = mat2np(inputs');
    inputs_now = inputs_now.reshape(int32(-1),int32(4));
    inputs_now = scaler_x.transform(inputs_now);
    
    outputs_now = NN.predict(inputs_now);
    outputs_now = scaler_y.inverse_transform(outputs_now);
    outputs_now = np2mat(outputs_now);
    outputs(:,1) = -outputs_now(:,1)/1000;
    outputs(:,2) = outputs_now(:,2)-param;
    %outputs(:,3) = -outputs_now(2)+48;
end

