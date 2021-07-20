function [fx,x,numeval]=dhc(fobj,x,initsize,thres,budget,x_L,x_U,weight,c_L,c_U,iterprint,tolc,varargin)

% n_out=nargout(fobj);

if ~isempty(c_L) | ~isempty(c_U)
    n_out=2;
else
    n_out=1;
end

NDIM=size(x,2);
THRESHOLD=thres;
INIT_SIZE=initsize;

numeval=0;

v=zeros(1,NDIM);        %Inicializar los vi con ceros
u=v;
vi=-1;
vvec=1;
vr=-INIT_SIZE;
xreal=x.*(x_U-x_L)+x_L;

if n_out>1
    [fx,cx]=feval(fobj,xreal,varargin{:});
    penalty=ssm_penalty_function(xreal,cx,c_L,c_U,tolc);
    fx=fx+weight*penalty;
else
    [fx]=feval(fobj,xreal,varargin{:});
    cx=0;
end

numeval=numeval+1;
fxv=1e30;

%JRB
nnn=0;

while (abs(vr)>=THRESHOLD)
    if (abs(vr)<2*THRESHOLD)
        maxiter=2*NDIM;
    else
        maxiter=2;
    end
    iter=0;
    while (fxv>=fx) & (iter<maxiter)
        if iter==0
            xv=x;
        else
            xv(vi+1)=xv(vi+1)-vr;               %OJO CON ESTO QUE NO LO TENGO CLARO
        end
        if (vvec)
            vvec=0;
        end
        vr=-vr;
        if (vr>0)
            vi=mod((vi+1),NDIM);
        end

        xv(vi+1)=xv(vi+1)+vr;               %OJO TB CON ESTO

        %Meter en bounds
        aaa=find(xv<0);
        bbb=find(xv>1);

        xv(aaa)=0;
        xv(bbb)=1;

        xvreal=xv.*(x_U-x_L)+x_L;

        if n_out>1
            [fxv,cxv]=feval(fobj,xvreal,varargin{:});
            penaltyv=ssm_penalty_function(xreal,cxv,c_L,c_U,tolc);
            fxv=fxv+weight*penaltyv;
        else
            [fxv]=feval(fobj,xvreal,varargin{:});
            cxv=0;
        end

        pen2=0;
        aaa=find(xvreal<x_L);
        bbb=find(xvreal>x_U);

        if ~isempty(aaa)
            pen2=pen2+sum((x_L(aaa)-xvreal(aaa)));
        end

        if ~isempty(bbb)
            pen2=pen2+sum((xvreal(bbb)-x_U(bbb)));
        end

        fxv=fxv+weight*pen2;
        
        numeval=numeval+1;
        iter=iter+1;
        
        if numeval>=budget
            vr=thres/10;
            break
        end

    end

    if (fxv>=fx) | isnan(fxv)
        fxv=1e30;
        vr=vr/2;
    else
        fx=fxv;
        x=xv;
        %JRB
        if iterprint
%             nnn=nnn+1;
%             plot(nnn,fx,'o');hold on;drawnow;
            fprintf('%s %i %s %g \n','Nevals:', numeval, 'Fobj:',fx)
        end
        %JRB

        if (iter==0)
            if (vvec)
                u=u+v;
                v=v*2;
                xv=xv+v;
                vr=vr*2;
            else
                u(vi+1)=u(vi+1)+vr;         %OJOOOOOO!!!!!
                vr=vr*2;
                xv(vi+1)=xv(vi+1)+vr;       %OJOOOOOOOOOO!!!!
            end
            %Meter en bounds
            aaa=find(xv<0);
            bbb=find(xv>1);

            xv(aaa)=0;
            xv(bbb)=1;

            xvreal=xv.*(x_U-x_L)+x_L;
            
            if n_out>1
                [fxv,cxv]=feval(fobj,xvreal,varargin{:});
                penaltyv=ssm_penalty_function(xreal,cxv,c_L,c_U,tolc);
                fxv=fxv+weight*penaltyv;
            else
                [fxv]=feval(fobj,xvreal,varargin{:});
                cxv=0;
            end

            pen2=0;
            aaa=find(xvreal<x_L);
            bbb=find(xvreal>x_U);

            if ~isempty(aaa)
                pen2=pen2+sum((x_L(aaa)-xvreal(aaa)));
            end

            if ~isempty(bbb)
                pen2=pen2+sum((xvreal(bbb)-x_U(bbb)));
            end

            fxv=fxv+weight*pen2;

            
            numeval=numeval+1;

            if numeval>=budget
                vr=thres/10;
                break
            end
        else
            xv=xv+u;
            xv(vi+1)=xv(vi+1)+vr;            %OJOOOOOOOOO!!!!

            %Meter en bounds
            aaa=find(xv<0);
            bbb=find(xv>1);

            xv(aaa)=0;
            xv(bbb)=1;

            xvreal=xv.*(x_U-x_L)+x_L;
            
            if n_out>1
                [fxv,cxv]=feval(fobj,xvreal,varargin{:});
                penaltyv=ssm_penalty_function(xreal,cxv,c_L,c_U,tolc);
                fxv=fxv+weight*penaltyv;
            else
                [fxv]=feval(fobj,xvreal,varargin{:});
                cxv=0;
            end

            pen2=0;
            aaa=find(xvreal<x_L);
            bbb=find(xvreal>x_U);

            if ~isempty(aaa)
                pen2=pen2+sum((x_L(aaa)-xvreal(aaa)));
            end

            if ~isempty(bbb)
                pen2=pen2+sum((xvreal(bbb)-x_U(bbb)));
            end

            fxv=fxv+weight*pen2;

            
            numeval=numeval+1;

            if numeval>=budget
                vr=thres/10;
                break
            end


            if (fxv>=fx) | isnan(fxv)
                u=zeros(1,NDIM);
                xv=x;
                u(vi+1)=vr;                  %OJOOOOOOOOO!!!!
                vr=vr*2;
                xv(vi+1)=xv(vi+1)+vr;                     %OJOOOOOOOOO!!!!

                %Meter en bounds
                aaa=find(xv<0);
                bbb=find(xv>1);

                xv(aaa)=0;
                xv(bbb)=1;

                xvreal=xv.*(x_U-x_L)+x_L;

                if n_out>1
                    [fxv,cxv]=feval(fobj,xvreal,varargin{:});
                    penaltyv=ssm_penalty_function(xreal,cxv,c_L,c_U,tolc);
                    fxv=fxv+weight*penaltyv;
                else
                    [fxv]=feval(fobj,xvreal,varargin{:});
                    cxv=0;
                end

                pen2=0;
                aaa=find(xvreal<x_L);
                bbb=find(xvreal>x_U);

                if ~isempty(aaa)
                    pen2=pen2+sum((x_L(aaa)-xvreal(aaa)));
                end

                if ~isempty(bbb)
                    pen2=pen2+sum((xvreal(bbb)-x_U(bbb)));
                end

                fxv=fxv+weight*pen2;
                numeval=numeval+1;
                
                if numeval>=budget
                    vr=thres/10;
                    break
                end

            else
                x=xv;
                fx=fxv;
                u(vi+1)=u(vi+1)+vr;         %OJOOOOOOOOOO!
                v=2*u;
                vvec=1;
                xv=xv+v;

                %Meter en bounds
                aaa=find(xv<0);
                bbb=find(xv>1);

                xv(aaa)=0;
                xv(bbb)=1;

                xvreal=xv.*(x_U-x_L)+x_L;
                
                if n_out>1
                    [fxv,cxv]=feval(fobj,xvreal,varargin{:});
                    penaltyv=ssm_penalty_function(xreal,cxv,c_L,c_U,tolc);
                    fxv=fxv+weight*penaltyv;
                else
                    [fxv]=feval(fobj,xvreal,varargin{:});
                    cxv=0;
                end

                pen2=0;
                aaa=find(xvreal<x_L);
                bbb=find(xvreal>x_U);

                if ~isempty(aaa)
                    pen2=pen2+sum((x_L(aaa)-xvreal(aaa)));
                end

                if ~isempty(bbb)
                    pen2=pen2+sum((xvreal(bbb)-x_U(bbb)));
                end

                fxv=fxv+weight*pen2;
                
                numeval=numeval+1;
                if numeval>=budget
                    vr=thres/10;
                    break
                end


                vr=0;
                vr=sum(v.^2);
                vr=sqrt(vr);
            end
        end
    end
end

if fxv<fx & ~isnan(fxv)
    fx=fxv;
    x=xv;
end

x=x.*(x_U-x_L)+x_L;