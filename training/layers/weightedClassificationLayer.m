classdef weightedClassificationLayer  < nnet.layer.ClassificationLayer
        
    properties
        % (Optional) Layer properties.
        ClassWeights
        % Layer properties go here.
    end
 
    methods
        function layer = weightedClassificationLayer (name,classWeights)           
            % (Optional) Create a myClassificationLayer.
            layer.ClassWeights = classWeights;
                        % Set layer name
            if nargin == 2
                layer.Name = name;
            end

            % Set layer description
            layer.Description = 'Weighted cross entropy';
            % Layer constructor function goes here.
        end

        function loss = forwardLoss(layer, Y, T)
            N = size(Y,4);
            Y = squeeze(Y);
            T = squeeze(T);
            W = layer.ClassWeights;
    
            loss = -sum(W*(T.*log(Y)))/N;
        end
        
    end
end