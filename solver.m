       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 INDIAN INSTITUTE OF TECHNOLOGY GUWAHATI                 %
%                  DEPARTMENT OF MECHANICAL ENGINEERING                   %
%                                                                         %
%                          2022-23 2ND SEMESTER                           %
%                                                                         %
%               ME 682 - NONLINEAR FINITE ELEMENT METHODS                 %
%                                                                         %
%                                                                         %
% Code initially developed by: Sachin Singh Gautam                        %
%                                                                         %
%                                                                         %
% Project 1: Due date 31.03.2023, Friday, 5 PM                            %
%                                                                         %
% The code is written for solving a finding the displacement, strains and % 
% stresses for a cantilever beam subjected to point load as shown below   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
I2 = eye(2) ; 
ksp = zeros(nksz,1) ; ndoel = ndoelo ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%                    Bulk Element Loop                                   %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kel_dummy = zeros(ndoel,ndoel,nel) ; rel_dummy = zeros(ndoel,nel) ;
jtermu = zeros(nel,1) ;
parfor i = 1:nel   % Run the loop over elements

   kel = zeros(ndoel,ndoel) ;% rel = zeros(ndoel,1) ; 

   [kel, rel, ie, jtermue] = get_element_stiffness_right_side_vector(...
       2*CON(i,:),ndoel,Xn(CON(i,:),:),xn(CON(i,:),:),U,ngpv,xigv,I2,D) ; 

    kel_dummy(:,:,i) = kel ; 
    rel_dummy(:,i) = rel ; % This is external laod vector due to distributed load.
    jtermu(i,1) = jtermue ; % Store error flag
end


face_con = [1 2 ; 2 3 ; 3 4 ; 4 1]; % fixed for every element
ndoele_surf = ny; % no. of dofs of the load surface
 ele_load_face = zeros(ndoele_surf,2) ;
for i=1:ndoele_surf
    ele_load_face(i,:) = [i*nx 2] ;
end

ele_load_face;
%ele_load_face = [6 3;7 3;8 3;9 3;10 3]; % first column = element no. ; second column = face no.
% nele_surface_load = nx; % no of elements having surface load
nele_surface_load = ny;
% t = [0 ; -1] ;
% t = [3000 ; 0] ;
    F_ext = zeros(ndof,1) ; 
    for i = 1:nele_surface_load
    element_no = ele_load_face(i,1) ;
    eleface = ele_load_face(i,2);
    local_con_face = face_con(eleface,:);
    global_con_face = CON(element_no,local_con_face);
    ig_1 = 2*global_con_face ;
    ie_global = [ig_1(1)-1 ig_1(1) ig_1(2)-1 ig_1(2)];
    ele_force_vector = zeros(4,1);
    if i <= nx/2
        t = [10000 ; 0] ;
    elseif i > nx/2
        t = [0 ; 15000] ;
    end
    for gp = 1:ngps
        xi = xigs(gp,1) ; wg = xigs(gp,2) ;
        N1 = (1-xi)/2 ; N2 = (1+xi)/2 ;
        N = [N1*I2 N2*I2] ;
        Jac = 0.5*((Xn(global_con_face(1),2)-Xn(global_con_face(2),2))^2 + (Xn(global_con_face(1),1)-Xn(global_con_face(2),1))^2)^0.5 ;
        ele_force_vector = ele_force_vector + N'*t*Jac*wg ;
    end
    F_ext(ie_global) = F_ext(ie_global) + ele_force_vector ;
     %rel_dummy(:,i+5) =  ele_force_vector ;
    end


% for i = 1:nel   
%    con = 2*CON(i,:) ;
%    ie = [con(1)-1 con(1) con(2)-1 con(2) con(3)-1 con(3) con(4)-1 con(4)] ;
%    F_ext(ie,1) = F_ext(ie,1) + rel_dummy(:,i) ;
% end

for i = 1:nel
  for j = 1:ndoel
   ksp(cspa_global(j,i):csp_global(j,i))  = kel_dummy(:,j,i) ;
  end
end

clear kel_dummy ;
clear rel_dummy ;

% Check if there was any error in the element loop then kindly throw the
% error message
if jterm == 0
    if sum(jtermu(:,1)) > 0
        jterm = 1 ;
        fprintf('\n\n Error in computation in element loop ..... exiting the newton iteration loop ...') ;
%             fprintf(fipf,'\n\n Error in computation in element loop ..... exiting the newton iteration loop ...') ; 
%             break ;        
    end
end

% Sparse Stiffness Matrix
K = sparse(isp,jsp,ksp) ;

clear ksp ;

% % Apply external point load 
% F_ext(dof_point_load,1) = -F_tip ;

% Eliminate Constraints in K and R due to BC and Condensation
Kr = K(ir,ir) ;
fr = F_ext(ir,1) ;

clear K ; clear f ;
   
% Get the sparse format K
Kr     = sparse(Kr) ;
  
% Solve for incremental displacement
DUr    =  Kr \ fr ;
   
clear Kr ;
   
DU     = zeros(ndof,1) ; % DU is the iterative displacement
DU(ir) = DUr ;  
U      = U + DU ; % Add incremental displacement to  previous displacement 
                  %to get the current total displacement
   
 % Current configuration
Ux = [ U(i1) U(i2) ] ; % Rearragne the displacement in the same form as coordinate array
xn = Xn + Ux(:,1:2) ;