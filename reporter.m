function []=reporter(x,y,n)
% test sample
%     x=([6 7 8 5 3 2 1 2 5 8 ; 9 7 0 8 5 6 4 3 2 10])';
%     y=([9 7 6 4 3 2 5 6 8 1])';
%     n=2;
    
    disp('-------- Multi Dimentional Multi Order Regression Solver --------')
    disp('-------------------------- Full Report --------------------------')
    
    sol=MainFunction(x,y,n);
    % report contents
    InputMatrix=[x,y]
    Order= n
    GeneratedMatrix = sol.nx
    Equation = sol.finalForm
    Properties = sol.props
end
