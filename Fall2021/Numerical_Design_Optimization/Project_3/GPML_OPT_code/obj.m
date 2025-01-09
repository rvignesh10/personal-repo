%% Objective function calculation

function f = obj(x)
    % x - design variable
    f = subObj(x);
    
    function fun = subObj(z)
        
        % set the squared exponential covariance function
        covfunc = {@covMaterniso,1};
        % trained hyper parameters - read from Workspace
        hyp_t = evalin('base','hyp');
        
        % set the likelihood function to Gaussian
        likfunc = @likGauss;
        
        % read training data 
        % x - training points; y - training obj function
        x_t = evalin('base','x');
        y_t = evalin('base','y');
       
        fun = gp(hyp_t, @infExact, [], covfunc, likfunc, x_t, y_t, z');
        
        % negate fun value to find the optimal point that generates
        % maximizer
        fun = -1 * fun;
    end

end