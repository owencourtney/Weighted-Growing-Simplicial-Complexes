function[sat]=check_saturation(nodes_nj,nodes_sel_tri)
% INPUTS:
% - nodes_nj, vector of the labels of the three nodes in the proposed
% triangle nj.
% - nodes_sel_tri, matrix of the labels of the nodes in the triangles
% selected so far.
% OUTPUT:
% - sat=0 if each link of the proposed triangle appears at most once in the
% in the triangles selected so far, sat=1 otherwise.
sat=0;
for i1=1:3
    for i2=i1+1:3
        counter=0;
        for n=1:length(nodes_sel_tri(:,1))
            if any(nodes_nj(i1)==nodes_sel_tri(n,:))&&any(nodes_nj(i2)==nodes_sel_tri(n,:))
                counter=counter+1;
            end
        end
        if counter>=2
            sat=1;
        end
    end
end

