classdef labelSmoothClassificationLayer  < nnet.layer.ClassificationLayer
        
    properties
        % (Optional) Layer properties.
        scale
        numClasses
        % Layer properties go here.
    end
 
    methods
        function layer = labelSmoothClassificationLayer (name,numClasses,scale)           
            % (Optional) Create a myClassificationLayer.
            layer.scale = scale;
            layer.numClasses = numClasses;
            % Set layer name
            layer.Name = name;
            if nargin == 2
                layer.Name = name;
            end

            % Set layer description
            layer.Description = 'labelSmooth cross entropy';
            % Layer constructor function goes here.
        end

        function loss = forwardLoss(layer, Y, T)
            N = size(Y,4);
            Y = squeeze(Y);
            T = squeeze(T);
            nC = layer.numClasses;
            scl = layer.scale;
            T = T * (1 - scl) + scl / nC;
            W = [1 1];
            loss = -sum(W*(T.*log(Y)))/N;
        end
        
    end
end