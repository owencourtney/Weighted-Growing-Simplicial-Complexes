%++++++++ Weighted Growing Simplicial Complexes in dimension 3++++++++++++++++++++++++++++++++++++++++++++
% 
% This code can be redistributed and/or modified
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%  
% This program is distributed ny the authors in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%  
% If you use this code please cite the following two papers:
%
%[1] O. T. Courtney and G. Bianconi, arXiv preprint, arXiv:1703.01187 (2017)
%[2] G. Bianconi, C. Rahmede  'Network geometry with flavor: From complexity to quantum geometry.' 
%	   Physical Review E 93 032315 (2016).
%
%
% (c) Owen.T.Courtney (email: o.t.o.courtney@qmul.ac.uk ) 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


function [a,tri_n,tet_n,tet_t,k_30,k_31,k_32,s_30,s_31,s_32,w,l_combined,t_combined] = WGSC_dim3(N,s,m,mp,h,figure)

% Code that generates Weighted Growing Simplicial Complexes (WGSCs) in dimension d=3 and flavours s=-1,0,1.

% Below we list the output variables and their meanings. We note that the
% full WGSC may be reconstructed from the variables tet_n and w.

% Outputs:
% a adjacency matrix.
% tri_n matrix of nodes (columns) of each triangle (rows).
% tet_n matrix of nodes (columns) of each tetrahedron (rows).
% tet_t matrix of triangles (columns) of each tetrahedron (rows).
% k_30 vector of generalized degrees k_{3,0} of the nodes.
% k_31 matrix of generalized degrees k_{3,1} of the links.  
% k_32 vector of generalized degrees k_{3,2} of the triangles.  
% s_30 vector of generalized degrees k_{3,0} of the nodes.
% s_31 matrix of generalized degrees k_{3,1} of the links.  
% s_32 vector of generalized degrees k_{3,2} of the triangles. 
% w vector of weights of the tetrahedra.
% l_combined combined link data - 1 row per link, 1st column node 1, 2nd
% column node 2, 3rd column generalized degree of link, 4th column
% generalized strength of link.
% t_combined combined triangle data - 1 row per triangle, 1st column node
% 1, 2nd column node 2, 3rd column node 3, 4th column generalized degree of
% triangle, 5th column generalized strength of triangle.

% This code uses: 
% N maximal number of nodes in the NGF
% Flavour of the WGSC s=-1,0,1
% m number of tetrahedra added at each timestep
% mp number of tetrahedra reinforced at each timestep
% h number of nodes in the 'early growth' phase
% figure=1 will save the weighted edge list of the skeleton-network with 
% link-weights given by the generalized strengths of the links in the SC.
% File name:
% "WGSC_edgelist_N%d_m%d_mp%d_s%d.csv"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
k_31=sparse(N,N);
s_31=sparse(N,N);
s_30=zeros(N,1);
k_30=zeros(N,1);
w=zeros(1,h-3 + m*(N-h));
nt=0; % Initialization number of triangles   
nt3=1; %Initialization number of tetrahedra

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Initial condition: at time t=1 a single tedrahedron (1,2,3,4)
w(nt3)=1;
for i1=1:4,
    tet_n(nt3,i1)=i1;
    k_30(i1)=1;
    s_30(i1)=1;
    for i2=(i1+1):4, 
        k_31(i1,i2)=1;
        k_31(i2,i1)=1; 
        s_31(i1,i2)=1;
        s_31(i2,i1)=1;
       for i3=(i2+1):4,           
           nt=nt+1;
           tet_t(nt3,nt)=nt;
           tri_n(nt,1)=i1;
           tri_n(nt,2)=i2;
           tri_n(nt,3)=i3;
           a_occ(nt)=1;
           k_32(nt)=1;
           s_32(nt)=1;
       end
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Early growth: nodes are added sequentially with each node contributing 1
% tetrahedra with weight 1, untill the simplicial complex contains h nodes.
% This early growth is equivalent to setting m=1, mp=0. Variable h should
% be chosen such that there are an adequate number of triangles and
% tetrahedra in the early complex for the 'growth' and 'reinforcement' 
% processes to function properly. This depends on the choices of m and mp.

for in=5:h,
    % select an existing triangle for the tetrahedra to attach to.
    [I,J,V]=find(a_occ);
    norm=sum(V);
    x=rand(1)*norm;
    for nj1=1:numel(V),
            x=x-V(nj1);
         if x<0,
             nj=J(nj1); 
             break; % select triangle with label nj
         end
    end
    
    % update number of tetrahedra
    nt3=nt3+1;
    
    % new node has label in and has initial strength and degree 1:
    s_30(in)=1;
    k_30(in)=1;
    
    % update occupancy, degree and strength of triangle nj:
    a_occ(nj)=a_occ(nj)+s;
    k_32(nj)=k_32(nj)+1;
    s_32(nj)=s_32(nj)+1;
    
    % retrieve labels of nodes in triangle nj:
    l(1)=tri_n(nj,1);
    l(2)=tri_n(nj,2);
    l(3)=tri_n(nj,3);
    
    % new tetrahedron is called nt3. Record labels of its nodes:
    for i1=1:3,
        tet_n(nt3,i1)=l(i1);
    end
    tet_n(nt3,4)=in;
    
    % new tetrahedron has initial weight 1:
    w(nt3)=1;
    
    % add new links from the new node to the existing triangle:
    for n=1:3,
        k_31(in,l(n))=1;
        k_31(l(n),in)=1;  
        s_31(in,l(n))=1;
        s_31(l(n),in)=1;  
    end
    
    % update degrees and strengths of nodes and links in existing triangle:
    for n1=1:3,
        s_30(l(n1))=s_30(l(n1))+1;
        k_30(l(n1))=k_30(l(n1))+1;
        for n2=n1+1:3,
            k_31(l(n1),l(n2))=k_31(l(n1),l(n2))+1;
            k_31(l(n2),l(n1))=k_31(l(n2),l(n1))+1;
            s_31(l(n1),l(n2))=s_31(l(n1),l(n2))+1;
            s_31(l(n2),l(n1))=s_31(l(n2),l(n1))+1;
        end
    end
    
    % add new triangles: update adjacencies, occupancy, degrees and
    % strengths
    tet_t(nt3,1)=nj;
    for n1=1:3,
        for n2=n1+1:3,           
            nt=nt+1;
            tet_t(nt3,n1+n2-1)=nt;
            tri_n(nt,1)=l(n1);
            tri_n(nt,2)=l(n2);
            tri_n(nt,3)=in;   
            a_occ(nt)=1;
            k_32(nt)=1;
            s_32(nt)=1;
        end
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main growth phase: simplicial complex now grows as specified in the
% paper. New node introduced at each time step, m tetrahedra created, mp
% tetrahedra reinforced.

for in=h+1:N,
    % Select m distinct triangles:
    nj_selected=zeros(size(a_occ));
    nj_rejected=zeros(size(a_occ));
    ntpause=nt;
    for mn=1:m,
        sat=1;
        while sat==1,
            [I,J,V]=find(a_occ.*(1-nj_selected).*(1-nj_rejected));
            % check there are enough available triangles to continue:
            if numel(V)<m-mn+1
                error('Insufficient available faces')
            end
            norm=sum(V);
            x=rand(1)*norm;
            for nj1=1:numel(V),
                x=x-V(nj1);
                if x<0,
                    nj=J(nj1);
                    break; % select triangle with label nj
                end
            end
            
            % if s=-1, check if proposed tetrahedra would create
            % 'oversaturated' triangles with triangles already selected:
            if s==-1
                nodes_nj=[tri_n(nj,1),tri_n(nj,2),tri_n(nj,3)];
                if mn==1
                    sat=0;
                else
                    [sat]=check_saturation(nodes_nj,nodes_sel_tri);
                    % sat=0 if no oversaturated triangles, sat=1 otherwise
                end
                if sat==1
                    % proposed tetrahedra would create illegal
                    % 'oversaturated' triangles. Reject nj and select from
                    % remaining tetrahedra:
                    nj_rejected(1,nj)=1; 
                end
            else
                % proposed tetrahedra is legal. Exit while loop.
                sat=0;
            end
        end
        
        if s==-1
            % create list of nodes of triangles selected so far:
            if mn==1
                nodes_sel_tri=nodes_nj;
            else
                nodes_sel_tri=[nodes_sel_tri;nodes_nj];
            end
        end
        % select triangle nj:
        nj_selected(1,nj)=1;
    end

    % create list of labels of selected triangles
    [~,selected_triangles]=find(nj_selected);

 
    % new node has label in and has initial strength and degree m:
    k_30(in)=m;
    s_30(in)=m;

    % add the m tetrahedra:
    for mn=1:m,
        % retrieve label of a selected exisiting triangle:
        nj=selected_triangles(mn);
        
        % update number of tetrahedra
        nt3=nt3+1;
        
        % update occupancy, degree and strength of triangle nj:
        a_occ(nj)=a_occ(nj)+s;
        k_32(nj)=k_32(nj)+1;
        s_32(nj)=s_32(nj)+1;
        
        % retrieve labels of nodes of the selected triangle
        l(1)=tri_n(nj,1);
        l(2)=tri_n(nj,2);
        l(3)=tri_n(nj,3);
        
        % new tetrahedron is called nt3. Record labels of its nodes:
        for i1=1:3,
            tet_n(nt3,i1)=l(i1);
        end
        tet_n(nt3,4)=in;
        
        % new tetrahedron has initial weight 1:
        w(nt3)=1;
        
        for n=1:3, % Add new links
            k_31(in,l(n))=k_31(in,l(n))+1;
            k_31(l(n),in)=k_31(l(n),in)+1;
            s_31(in,l(n))=k_31(in,l(n));
            s_31(l(n),in)=k_31(l(n),in);
        end
    
        % update degrees and strengths of nodes and links in existing triangle:
        for n1=1:3,
            s_30(l(n1))=s_30(l(n1))+1;
            k_30(l(n1))=k_30(l(n1))+1;
            for n2=n1+1:3,
                k_31(l(n1),l(n2))=k_31(l(n1),l(n2))+1;
                k_31(l(n2),l(n1))=k_31(l(n2),l(n1))+1;
                s_31(l(n1),l(n2))=s_31(l(n1),l(n2))+1;
                s_31(l(n2),l(n1))=s_31(l(n2),l(n1))+1;
            end
        end
        
        % Update: no. triangles (nt), tet-tri adjacency (As), triangle node-lists (tri_n), triangle energies (at), triangle occupancy expresions (a_occ), gen degrees of triangles (k_32)
        tet_t(nt3,1)=nj;
        for n1=1:3, 
            for n2=n1+1:3, 
                % check if the triangle formed by new node in and exisiting
                % link (l(n),(l(n2)) has already been added as part of a
                % previous tetrahedra:
                red=0;
                for ptri=(ntpause+1):nt,
                    if l(n1)==tri_n(ptri,1)&&l(n2)==tri_n(ptri,2),
                        red=1;
                        ftri=ptri; % triangle has already been added and has label ftri
                        break;
                    end
                end
                if red==0,
                    % Triangle has not be added previously. Add triangle:
                    nt=nt+1;
                    tet_t(nt3,n1+n2-1)=nt;
                    tri_n(nt,1)=l(n1);
                    tri_n(nt,2)=l(n2);
                    tri_n(nt,3)=in;   
                    a_occ(nt)=1;  
                    k_32(nt)=1;
                    s_32(nt)=1;
                else
                    % Triangle has been added previously. Update triangle:
                    tet_t(nt3,n1+n2-1)=ftri;
                    a_occ(ftri)=a_occ(ftri)+s;
                    k_32(ftri)=k_32(ftri)+1;
                    s_32(ftri)=s_32(ftri)+1;
                end
            end
        end  
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reinforce the tetrahedra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Choose the number of tetrahedra that will be reinforced (nmp):
    nmp=mp;
    if nt3<nmp,
            nmp=nt3;
    end
    
    % Randomly select nmp tetrahdera and reinforce them:
    pnj_selected=zeros(size(w));
    for itn=1:nmp,
    	[~,J,V]=find(w.*(1-pnj_selected));
    	norm=sum(V);
    	x=norm*rand(1);
        for ia=1:numel(V),
        	x=x-V(ia);
        	if (x<0)
                pnj=ia;
            	break; % select tetrahedron pnj
            end
        end
        pnj_selected(J(pnj))=1;
    end
            
	[~,selected_tetrahedra]=find(pnj_selected);

	for itn=1:nmp,
        % Retrieve label of a selected tetrahedron:
        inn=selected_tetrahedra(itn);

        % Update weight of tetrahedron inn:
        w(inn)=w(inn)+1;
            
        % Update strengths of triangles of inn:
        for i=1:4,
            nj=tet_t(inn,i);
            s_32(nj)=s_32(nj)+1;
        end
        
        % Update strengths of the nodes and links of inn:
        for i1=1:4,
            s_30(tet_n(inn,i1))=s_30(tet_n(inn,i1))+1;
            for i2=i1+1:4,
                s_31(tet_n(inn,i1),tet_n(inn,i2))=s_31(tet_n(inn,i1),tet_n(inn,i2))+1;
                s_31(tet_n(inn,i2),tet_n(inn,i1))=s_31(tet_n(inn,i2),tet_n(inn,i1))+1;
            end
        end

	end
            
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generalized degrees and strengths
%k_30=sum(k_31>0);
%k_30=k_30-2;
[I2,J2,A2]=find(tril(k_31));
[I3,J3,A3]=find(tril(s_31));
l_combined=[I2,J2,A2,A3];
a=k_31>0;
for ntri=1:nt,
    t_combined(ntri,:)=[ntri,tri_n(ntri,1),tri_n(ntri,2),tri_n(ntri,3),k_32(ntri),s_32(ntri)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print skeleton network file
if figure==1
    filename=sprintf('WGSC_edgelist_N%d_m%d_mp%d_s%d.csv',N,m,mp,s);
    fid=fopen(filename,'w');
    fprintf(fid, '%s\t%s\t%s\t%s\n','Source','Target','Weight','Type');
    for it=1:length(k_31),
        fprintf(fid, '%d\t%d\t%d\t%s\n', l_combined(it,1), l_combined(it,2),l_combined(it,4),'Undirected');
    end
    fclose(fid);
end



