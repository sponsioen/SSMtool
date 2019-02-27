function [Achar] = sym2char(A,ndof)

    
    symbolic = strcmp(class(A),'sym');

    if symbolic

    Mcell = sym2cell(A);
    Mchar = '[';
    
        for i = 1:ndof(1)
           for j = 1:ndof(2)
%               Mchar = strcat(Mchar,sprintf('%s',Mcell{i,j}));
                Mchar = strcat(Mchar,char(Mcell{i,j}));
              if j == ndof(2) && i ~= ndof(1)
                Mchar = strcat(Mchar,';');
              elseif  j == ndof(2) && i == ndof(1)
                Mchar = strcat(Mchar,']');
              else
                Mchar = strcat(Mchar,',');
              end
           end
        end

    else
        
    M = A;
    Mchar = '[';
    
        for i = 1:ndof(1)
           for j = 1:ndof(2)
              Mchar = strcat(Mchar,sprintf('%s',num2str(M(i,j))));

              if j == ndof(2) && i ~= ndof(1)
                Mchar = strcat(Mchar,';');
              elseif  j == ndof(2) && i == ndof(1)
                Mchar = strcat(Mchar,']');
              else
                Mchar = strcat(Mchar,',');
              end
           end
        end
        
        
    end
    Achar = Mchar; 
end 