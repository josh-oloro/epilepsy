Y_Pred_all = [];
window_size = 500;
window_skip = 1;

% layer_name = {'CA1'; 'Interlayer'; 'DG'}
layer_name = {'CA1'}

for l = 1:length(layer_name)

    wf = wave.(layer_name{l}).shf;
    t_max = length(wf);
%     t_max = 1500
    
    %Sliding window
    for wi = 1:window_skip:(t_max-window_size)
        wj = wi + window_size;
        window = wf(wi:wj-1);
        window_fsst = vertcat(window',real(fsst(window,10,kaiser(100,10))));
 
%         DataParts = cell(size(window_fsst,1), size(window_fsst,2), 2);
% 
%         for k = 1:length(window_fsst)
%             DataParts{k,:,1} = real(window_fsst{k,1});
%             %DataParts{k,:,2} = imag(XTrain{k,1});
% 
%             window_fsst{k,:} = DataParts{k,:,1};
%             window_fsst{k,:} = vertcat(window_fsst{k,:},DataParts{k,:,2});
%         end



        Y_Pred = classify(net,window_fsst);
        Y_Pred = repmat(Y_Pred,window_skip,1); % Verify
        Y_Pred_all = [Y_Pred_all,Y_Pred];

    end
end


figure;
fig1 = subplot(2,1,1)
plot(wave.CA1.shf)
fig2 = subplot(2,1,2)
plot(Y_Pred_all)
linkaxes([fig1,fig2],'x')
