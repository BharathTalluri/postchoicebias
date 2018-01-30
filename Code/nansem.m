function [y,n]=nansem(x,dim,num)

%now modified so that it takes a minimum number of rows parameter

%gives the nansem
%now also works for x is a matrix

if ~exist('num','var'),num=1;end

if size(x,1)==1&&~(exist('dim','var')), y=nanstd(x)/sqrt(nanlength(x)-1);
else
%elseif size(x,1)>1
    if ~exist('dim','var')||isnan(dim)||isempty(dim),dim=1;end
    foo=ones(size(x,1),size(x,2));
    foo(isnan(x))=0;
    n=nansum(foo,dim);
    yo=n>=num;
    y=nanstd(x,1,dim)./sqrt(n-1);
    if dim==1
    y(:,~yo)=NaN;
    else
    y(~yo,:)=NaN;    
    end
end

