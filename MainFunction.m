% Multi Order Multi Dimentional Regression Solver
% By Jamal Sa'd, Mechatronics Engineering, JU.
% for Numerical Engineering Methods Course
% Prof. Zaer Abu Hamour
% x -> input matrix, y -> output matrix, n-> regression order

function [ solution ] = MainFunction(x,y,n)

% % Sample 1D Data
% x=[1 2 3 4 5];
% y=[3 2 0.5 1.5 4];
% % Sample 2D Data
% x=[6 7 8 5 3 2 1 2 5 8 ; 9 7 0 8 5 6 4 3 2 10];
% y=[9 7 6 4 3 2 5 6 8 1];
%  y=y';
%  x=x';
%  n=3;


[O,F]=size(x);
% create combinations of all new generated features
Fmat = createCombinations(n,O,F);
% sort new features in an increasing order
Fmat = prettyMatrix(Fmat);
% generate strings for the new features
Fstr= createString(Fmat);
% create a new matrix for the new generated order-based features
nx = createNewX(x,Fmat);
% solve the regression formula using normal equation
coefs=((nx)'*(nx))\((nx)'*y);
% generate a string for the final formula
strco= strcat(num2str(coefs(:)),Fstr(:));

% text processing
strco=deblank(strco);
strco=reverse(deblank(reverse(strco)));
for i = 2:size(strco,1)
    q=char(strco(i));
    k=strco(i);
    if(q(1)~='-')
    k= strcat('+',k);
    end
    strco(i)=k;
end
fstrco=string(strjoin(strco));

% create symbolic equation for plotting
xnames='@(';
for i =1:F
    xnames= [xnames,'x',num2str(i)];
    if i ~= F
        xnames= [xnames,','];
    end
end
xnames=[xnames,')'];
finalEq=[xnames,fstrco];
newFstr=regexprep(Fstr,'[.*]','');



% return generated regression outputs as a struct object

solution.Fstr =newFstr;
solution.nx = nx;
solution.coefs = coefs;
solution.finalForm= fstrco;
solution.equation = str2func(strjoin(finalEq));


% Find properities of the new equation
propStr = findProperities(x,y,solution.equation);
solution.props= propStr;


% plot functions, 2D plane for 2D data, 3D plane for 3D data and multi
% scatter plot for higher order data

Draw();




function [propStr]= findProperities(x,y,eq)
    ymean = mean(y);
    xcells = num2cell(x);
    yestimated = zeros(size(xcells,1),1);
    for i=1:size(xcells,1)
        yestimated(i) = eq(xcells{i,:});
    end
    Sr = sum((y-yestimated).^2);
    St = sum((y-ymean).^2);
    R2 =(St-Sr)/St;
    propStr = strcat('Sr = ',num2str(Sr),', St = ',num2str(St),', R^2 = ',num2str(R2));
end
    
    
function []=Draw()
  
    if(F==1)
        X=x;
        Y=y;
        fx = min(X):(max(X)-min(X))/50:max(X)+0.2*max(X);
        yx = solution.equation(fx);

        plot(fx,yx,'LineWidth',3)
        grid on
        hold on
        plot(X,Y,'X','MarkerSize',10,'LineWidth',3) 
        
    end
    if (F == 2)
            plot3([1],[1],[1])
            hold on
            for i=1:1:O
                X=x';
                Y=y;
                X1=X(1,i);
                X2=X(2,i);
                t = solution.equation(X1,X2);
                if(Y(i) > t)
                    plot3(X(1,i),X(2,i),Y(i),'O','Color','r','MarkerSize',10,'MarkerFaceColor','r')
                    plot3(X(1,i),X(2,i),t:0.1:Y(i),'O','Color','r','MarkerSize',1,'MarkerFaceColor','r')
                else
                    plot3(X(1,i),X(2,i),Y(i),'O','Color','b','MarkerSize',10,'MarkerFaceColor','b')
                    plot3(X(1,i),X(2,i),Y(i):0.1:t,'O','Color','b','MarkerSize',1,'MarkerFaceColor','b')
                end

            end


            xlabel('X1')
            ylabel('X2')
            zlabel('Y')
            title('plot of Multiple Regression Equation')
            grid on
                    mn1=min(X(1,:));
                    mx1=max(X(1,:));
                    mn2=min(X(2,:));
                    mx2=max(X(2,:));
                    mnf=min(mn1,mn2);
                    mxf=max(mx1,mx2);
            [X1,X2] = meshgrid(mnf: (mxf-mnf)/50:mxf,mnf: (mxf-mnf)/50:mxf);
            Z = solution.equation(X1,X2);
            surf(X1,X2,Z,'FaceAlpha',0.4);
            hold off
        
    end
    
    if(F > 2)
       gplotmatrix(x,y,y);
      
    end
    

end

function [nx]= createNewX(x,Fmat)
[NF,OF]=size(Fmat);
O=size(x,1);
nx=zeros(O,NF);
for i=1:NF
    k=ones(O,1);
    for j=1:OF
        k=k.*( x(:,j).^(Fmat(i,j)) );
    end
    nx(:,i) =k;
end
end

function [Fstr] = createString(Fmat)
[r,c]=size(Fmat);

for i=1:r
    str='';
    for j=1:c
        if Fmat(i,j) >0
          str= strcat(str,'.*x',num2str(j),'.^',num2str(Fmat(i,j)));
        end
    end
k=strcat(str,'');
Fstr(i)=string(k);
Fstr(1)= '';
end
end

function [Fmat] = createCombinations(n,O,F)
Fmat=zeros(n^F,F);
for i=1:F
    Fmat(:,i)=floor(mod((1:n^F)/n^(i-1),n));
end
Fmat=[Fmat;eye(F)*n];
Fmat=Fmat(sum(Fmat,2)<=n,:);
end

% pretty Matrix function uses bubble sort algorithm to arrange features
% combinations from least to highest

function [Fmat] = prettyMatrix(Fmat)
[r,c]=size(Fmat);
Fmat(:,c+1)=zeros(r,1);
c=c+1;
for i=1:r
    Fmat(i,c)=sum(Fmat(i,[1:c]));
end
for j=1:r
    for i=1:r-1
        if (Fmat(i,c) > Fmat(i+1,c))
            Fmat(i,c);
            temp=Fmat(i+1,:);
            Fmat(i+1,:)=Fmat(i,:);
            Fmat(i,:)=temp;
        end
    end
end
Fmat(:,c) = [];
end

end